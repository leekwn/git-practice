#!/usr/bin/env python3
"""
1) Load aligned influenza A segment sequences (PB2, PB1, PA, HA, NP, NA, M, NS)
2) Keep only isolates that have all 8 segments
3) Cluster each segment into lineages (PB2-A, PB2-B, ...)
4) Define genotype for each isolate as the combination of 8 segment lineages
"""

from Bio import SeqIO
import os
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

# --- CONFIGURATION ---

ALIGNED_DIR = "aligned_segments"

# Expected aligned FASTA filenames
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
ALN_FILES = {seg: os.path.join(ALIGNED_DIR, f"{seg}_H6N1_world_2010_2025_aln2.fasta") for seg in SEGMENTS}

# --- HELPER FUNCTIONS ---

def get_isolate_id(rec_id):
    """
    Extract isolate name from a FASTA record ID.

    Strategy:
      1. Take the first token before whitespace.
      2. Then take the part before the first '|'.

    Examples:
      'A/chicken/Korea/123/2018|PB2|EPI123456' -> 'A/chicken/Korea/123/2018'
      'A/duck/Taiwan/01/2015 PB1 segment'      -> 'A/duck/Taiwan/01/2015'
    """
    s = rec_id
    #s = s.split()[0]      # before any whitespace
    s = s.split("|")[0]   # before any pipe
    return s
#new
def load_segment_sequences(seg, filepath):
    """
    Load aligned sequences for a given segment into a dict:
      { isolate_id: sequence_string }

    If an isolate appears multiple times in this segment,
    keep only the FIRST occurrence and ignore the others.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Aligned file for {seg} not found: {filepath}")

    seq_dict = {}
    dup_counts = {}

    for rec in SeqIO.parse(filepath, "fasta"):
        isolate = get_isolate_id(rec.id)
        seq_str = str(rec.seq).upper()

        if isolate in seq_dict:
            # already saw this isolate → skip this one
            dup_counts[isolate] = dup_counts.get(isolate, 1) + 1
            continue

        # first time we see this isolate → keep it
        seq_dict[isolate] = seq_str

    if dup_counts:
        total_dups = sum(c - 1 for c in dup_counts.values())
        print(
            f"{seg}: found {len(dup_counts)} isolates with duplicates "
            f"({total_dups} extra sequences skipped). "
            f"Examples: {list(dup_counts.keys())[:5]}"
        )

    print(f"{seg}: loaded {len(seq_dict)} unique isolates from {os.path.basename(filepath)}")
    return seq_dict


def hamming_distance(s1, s2):
    """
    Fraction of different positions between two aligned sequences.
    """
    assert len(s1) == len(s2), "Sequences must be aligned and same length."
    arr1 = np.frombuffer(s1.encode("ascii"), dtype="S1")
    arr2 = np.frombuffer(s2.encode("ascii"), dtype="S1")
    return (arr1 != arr2).mean()

def build_distance_matrix(seqs):
    """
    Build an NxN pairwise Hamming distance matrix for a list of sequences.
    """
    n = len(seqs)
    dist = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            d = hamming_distance(seqs[i], seqs[j])
            dist[i, j] = d
            dist[j, i] = d
    return dist

def cluster_segment(df_seg, max_distance=0.05):
    """
    Cluster sequences for a single segment into lineages.

    df_seg: DataFrame with columns ['isolate', 'segment', 'sequence'] for one segment.
    max_distance: rough target for maximum within-cluster distance (Hamming fraction).

    Returns: pd.Series of lineage labels, index = df_seg.index
    """
    seqs = df_seg["sequence"].tolist()
    isolates = df_seg["isolate"].tolist()
    n = len(seqs)

    if n == 1:
        seg = df_seg["segment"].iloc[0]
        return pd.Series([f"{seg}-A"], index=df_seg.index)

    # Build distance matrix
    dist_matrix = build_distance_matrix(seqs)

    # Estimate number of clusters from distance distribution (very rough heuristic)
    upper_tri = dist_matrix[np.triu_indices(n, k=1)]
    if len(upper_tri) == 0:
        n_clusters = 1
    else:
        max_d = upper_tri.max()
        # Avoid zero or negative
        n_clusters = max(1, int(round(max_d / max_distance)) or 1)

    seg = df_seg["segment"].iloc[0]
    print(f"{seg}: clustering {n} sequences into ~{n_clusters} clusters (max pairwise d ≈ {max_d:.3f})")
    
    model = AgglomerativeClustering(
        n_clusters=n_clusters,
        metric="precomputed",   # NEW name
        linkage="average"
    )
    #model = AgglomerativeClustering(
    #    n_clusters=n_clusters,
    #    affinity="precomputed",
    #    linkage="average"
    #)
    labels = model.fit_predict(dist_matrix)

    # Convert numeric labels 0,1,2,... into A,B,C,... for readability
    lineage_labels = []
    for lab in labels:
        lineage_labels.append(f"{seg}-{chr(ord('A') + lab)}")

    return pd.Series(lineage_labels, index=df_seg.index)

# --- MAIN PIPELINE ---

def main():
    # 1) Load all segments into dicts
    seg_to_seqs = {}
    for seg in SEGMENTS:
        seg_to_seqs[seg] = load_segment_sequences(seg, ALN_FILES[seg])

    # 2) Keep only isolates that have all 8 segments
    sets_of_isolates = [set(d.keys()) for d in seg_to_seqs.values()]
    common_isolates = set.intersection(*sets_of_isolates)
    print(f"\nIsolates with all 8 segments: {len(common_isolates)}")

    # 3) Build a tidy DataFrame: one row per (isolate, segment)
    rows = []
    for seg in SEGMENTS:
        for iso in common_isolates:
            seq = seg_to_seqs[seg].get(iso)
            if seq is None:
                # theoretically shouldn't happen because of the intersection
                continue
            rows.append({"isolate": iso, "segment": seg, "sequence": seq})

    seq_df = pd.DataFrame(rows)
    print(f"Total rows in seq_df (isolate x segment): {len(seq_df)}")

    # 4) Cluster each segment into lineages
    lineage_dfs = []
    for seg in SEGMENTS:
        df_seg = seq_df[seq_df["segment"] == seg].copy()
        df_seg["lineage"] = cluster_segment(df_seg, max_distance=0.01)
        lineage_dfs.append(df_seg[["isolate", "segment", "lineage"]])

    lineage_df = pd.concat(lineage_dfs, ignore_index=True)

    # 5) Pivot to get one row per isolate, columns = segment lineages
    geno_table = lineage_df.pivot_table(
        index="isolate",
        columns="segment",
        values="lineage",
        aggfunc="first"
    ).reset_index()

    # 6) Build genotype string (ordered by SEGMENTS list)
    def genotype_string(row):
        parts = []
        for seg in SEGMENTS:
            parts.append(row.get(seg, "NA"))
        return "/".join(parts)

    geno_table["genotype"] = geno_table.apply(genotype_string, axis=1)

    # 7) Save to CSV
    out_csv = "H6N1_genotypes.csv"
    geno_table.to_csv(out_csv, index=False)
    print(f"\nSaved genotype table to {out_csv}")
    print("Columns:", geno_table.columns.tolist())
    print("Example rows:")
    print(geno_table.head())

if __name__ == "__main__":
    main()
