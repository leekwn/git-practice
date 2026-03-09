#!/usr/bin/env python3
"""
Compute Korea x Japan seasonal pairwise HA identity counts >= 0.99 and plot heatmap.

Outputs:
 - counts_matrix.csv  (rows: Korea seasons, cols: Japan seasons)
 - pct_matrix.csv     (same, proportion of comparisons that >= 0.99)
 - heatmap_counts.png
 - heatmap_pct.png
"""

import argparse
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import re

# ---------- Helper utilities ----------
def parse_date(d):
    # robust date parser: try pandas
    if pd.isna(d):
        return None
    try:
        return pd.to_datetime(d, errors='coerce')
    except Exception:
        return pd.to_datetime(str(d), errors='coerce')

def season_label_from_date(dt):
    # season defined July 1 YEAR -> June 30 (YEAR+1)
    # e.g., date 2020-07-01..2021-06-30 -> "20/21"
    if dt is None or pd.isna(dt):
        return None
    year = dt.year
    # if month >= 9, season starts that year
    if dt.month >= 9:
        start = year
    else:
        start = year - 1
    end = start + 1
    return f"{str(start)[-2:]}/{str(end)[-2:]}"

def normalize_country(c):
    if pd.isna(c):
        return None
    s = str(c).strip().lower()
    if "korea" in s or s in ("kr","republic of korea","south korea","korea, south"):
        return "Korea"
    if "japan" in s or s in ("jp",):
        return "Japan"
    return None

def find_metadata_key_in_header(header, keys):
    # Try to match accession or virus-name tokens in the header
    hdr = header
    for k in keys:
        # require token match: word boundary or exact substring for EPI_ISL
        if re.search(r'(?<!\w){}\b'.format(re.escape(k)), hdr):
            return k
        # also try simple substring as fallback
        if k in hdr:
            return k
    return None

def pairwise_identity(seqA, seqB):
    # compute identity over positions where neither is gap (-)
    # seqA & seqB are aligned strings (same length)
    if len(seqA) != len(seqB):
        raise ValueError("Sequences must be same length for aligned identity")
    matches = 0
    compared = 0
    for a,b in zip(seqA, seqB):
        if a == '-' or b == '-':
            continue
        compared += 1
        if a == b:
            matches += 1
    if compared == 0:
        return 0.0
    return matches / compared

# ---------- Main pipeline ----------
def main(args):
    # 1) read metadata
    meta = pd.read_csv(args.meta, dtype=str)
    # Normalize column names (lowercase)
    meta_cols = {c.lower():c for c in meta.columns}
    # find accession/virus and date and country columns
    acc_col = None
    for candidate in ("isolate_id","accession"):#,"accession id","epi_isl","epi_isl_id","gisaid_epi_isl","strain","virus name","virus_name","virus"):
        if candidate in meta_cols:
            acc_col = meta_cols[candidate]
            break
    date_col = None
    for candidate in ("collection_date","date","collection date"):
        if candidate in meta_cols:
            date_col = meta_cols[candidate]
            break
    country_col = None
    for candidate in ("country","location","host_location","country_region"):
        if candidate in meta_cols:
            country_col = meta_cols[candidate]
            break

    if acc_col is None:
        raise SystemExit("Could not find accession/strain column in metadata. Please include a column named 'accession' or 'virus name' (case-insensitive).")

    if date_col is None:
        raise SystemExit("Could not find collection date column in metadata (collection_date).")

    if country_col is None:
        raise SystemExit("Could not find country column in metadata (country).")

    # parse dates & season & country
    meta['parsed_date'] = meta[date_col].apply(parse_date)
    meta['season'] = meta['parsed_date'].apply(season_label_from_date)
    meta['country_norm'] = meta[country_col].apply(normalize_country)

    # Build lookup mapping from accession/strain token to metadata row (allow duplicates by picking first)
    # Trim whitespace from keys
    key_to_row = {}
    for i,row in meta.iterrows():
        key = str(row[acc_col]).strip()
        if key and key not in key_to_row:
            key_to_row[key] = row

    print(f"Metadata loaded: {len(meta)} rows, keys for lookup: {len(key_to_row)}")
    print("Seasons present:", sorted(meta['season'].dropna().unique()))
    print("Countries present (normalized):", sorted(meta['country_norm'].dropna().unique()))

    # 2) load aligned sequences and attach metadata
    seq_records = list(SeqIO.parse(args.align, "fasta"))
    print(f"Aligned sequences loaded: {len(seq_records)}")
    # Build a mapping from sequence record -> (country, season)
    seq_info = {}   # header -> dict
    # For matching, try to match accession keys in header tokens
    keys = list(key_to_row.keys())

    for rec in seq_records:
        hdr = rec.description
        found_key = find_metadata_key_in_header(hdr, keys)
        if found_key:
            row = key_to_row[found_key]
            country = normalize_country(row[country_col])
            season = season_label_from_date(parse_date(row[date_col]))
            seq_info[rec.id] = {
                "header": hdr,
                "country": country,
                "season": season,
                "seq": str(rec.seq).upper()
            }
        else:
            # try looser matching: case-insensitive substring of virus name if separate column exists
            # Look for virus name column as fallback
            fallback_matched = False
            for candidate_col in [c for c in meta.columns if 'virus' in c.lower() or 'strain' in c.lower()]:
                val = str(meta.loc[meta.index, candidate_col]) if False else None
            # mark unmatched
            seq_info[rec.id] = {"header": hdr, "country": None, "season": None, "seq": str(rec.seq).upper()}
    # count assigned
    assigned = sum(1 for v in seq_info.values() if v['country'] in ("Korea","Japan") and v['season'] is not None)
    print(f"Assigned Korea/Japan + season to {assigned} / {len(seq_records)} sequences")

    # 3) group sequences by country and season
    groups = defaultdict(list)  # key: (country, season) -> list of seq strings
    for info in seq_info.values():
        c, s = info['country'], info['season']
        if c in ("Korea","Japan") and s is not None:
            groups[(c,s)].append(info['seq'])

    # build lists of Korea seasons and Japan seasons seen
    korea_seasons = sorted({s for (c,s) in groups.keys() if c=="Korea"})
    japan_seasons = sorted({s for (c,s) in groups.keys() if c=="Japan"})

    print("Korea seasons:", korea_seasons)
    print("Japan seasons:", japan_seasons)

    if not korea_seasons or not japan_seasons:
        raise SystemExit("No sequences found for Korea or Japan with season information. Check metadata matching and date parsing.")

    # 4) compute pairwise counts and proportions
    # initialize matrices
    counts = pd.DataFrame(0, index=korea_seasons, columns=japan_seasons, dtype=int)
    totals = pd.DataFrame(0, index=korea_seasons, columns=japan_seasons, dtype=int)

    threshold = args.threshold  # e.g., 0.99

    for ks in korea_seasons:
        seqs_k = groups.get(("Korea", ks), [])
        for js in japan_seasons:
            seqs_j = groups.get(("Japan", js), [])
            # total comparisons = len(seqs_k) * len(seqs_j)
            total_pairs = len(seqs_k) * len(seqs_j)
            totals.loc[ks, js] = total_pairs
            if total_pairs == 0:
                counts.loc[ks, js] = 0
                continue

            # compute pairwise identities
            above = 0
            # nested loops; optimize by early repr if huge - but we do brute force here
            for a in seqs_k:
                for b in seqs_j:
                    iden = pairwise_identity(a, b)
                    if iden >= threshold:
                        above += 1
            counts.loc[ks, js] = above
            print(f"{ks} (Korea) x {js} (Japan): total={total_pairs}, >={threshold*100:.2f}% id={above}")

    # compute percentage matrix
    pct = counts.copy().astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):
        pct = pct / totals.replace({0: np.nan})

    # 5) output CSVs
    counts.to_csv("counts_matrix.csv")
    pct.to_csv("pct_matrix.csv")

    # 6) plot heatmaps (counts and percentages)
    plt.figure(figsize=(max(6, len(japan_seasons)*0.8), max(6, len(korea_seasons)*0.6)))
    sns.heatmap(counts, annot=True, fmt="d", cmap="YlGnBu", cbar_kws={'label': 'count pairs >= threshold'})
    plt.title(f"Korea x Japan HA pairwise counts (identity >= {threshold*100:.1f}%)")
    plt.ylabel("Korea season")
    plt.xlabel("Japan season")
    plt.tight_layout()
    plt.savefig("heatmap_counts.png", dpi=300)
    plt.close()

    plt.figure(figsize=(max(6, len(japan_seasons)*0.8), max(6, len(korea_seasons)*0.6)))
    sns.heatmap(pct, annot=True, fmt=".3f", cmap="YlOrRd", cbar_kws={'label': 'fraction pairs >= threshold'}, vmin=0, vmax=1)
    plt.title(f"Korea x Japan HA pairwise fraction (identity >= {threshold*100:.1f}%)")
    plt.ylabel("Korea season")
    plt.xlabel("Japan season")
    plt.tight_layout()
    plt.savefig("heatmap_pct.png", dpi=300)
    plt.close()

    print("Done. Outputs:")
    print(" - counts_matrix.csv")
    print(" - pct_matrix.csv")
    print(" - heatmap_counts.png")
    print(" - heatmap_pct.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Korea-Japan seasonal HA pairwise identity counts and heatmaps")
    parser.add_argument("--align", required=True, help="Aligned HA fasta (MAFFT output)")
    parser.add_argument("--meta", required=True, help="Metadata CSV with accession/strain, collection_date, country")
    parser.add_argument("--threshold", type=float, default=0.99, help="Identity threshold (e.g., 0.99)")
    args = parser.parse_args()
    main(args)
