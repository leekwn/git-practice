"""
Agent-style pipeline to retrieve HA sequences for clade 2.3.4.4b HPAI (last 4 years),
extract the HA cleavage-site region, create a multiple-sequence alignment, and
produce summary tables of cleavage-site motifs.

Notes & caveats:
- GISAID holds the most comprehensive and up-to-date H5 datasets. Many sequences
  and clade annotations are present in GISAID but are restricted by its terms of use.
  This script provides a GISAID placeholder function: you'll need to download sequences
  yourself (or obtain API access) and place them in FASTA format, or provide credentials
  and modify the placeholder.
- NCBI GenBank contains many HA sequences and is accessible via Entrez — but clade
  metadata may be inconsistent. Filtering by collection date and manual curation
  is often necessary.
- MAFFT must be installed in your PATH for the alignment step.
- Python packages required: biopython, pandas, tqdm.

Outputs:
- data/raw_genbank_ha.fasta         -- raw HA nucleotide sequences fetched from GenBank
- data/gisaid_ha_placeholder.fasta  -- (optional) user-provided GISAID sequences
- data/combined_ha.fasta            -- combined nt sequences (deduplicated)
- data/ha_cleavage_regions.fasta    -- extracted nucleotide sequences around cleavage site
- results/ha_cleavage_alignment.fasta -- MAFFT amino-acid alignment of cleavage-region
- results/cleavage_motif_table.csv  -- table of accession, collection_date, motif, comments

Usage example:
    python agent_h5_clade_2.3.4.4b_pipeline.py --email you@example.com \
        --start 2021-09-14 --end 2025-09-14 --outdir ./results

"""

import os
import sys
import argparse
import subprocess
from datetime import datetime
from collections import defaultdict

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
from tqdm import tqdm

# --------------------------- Configuration ---------------------------------
# Default search parameters. Adjust as needed.
DEFAULT_START = "2021-09-14"  # four years prior to 2025-09-14 as requested
DEFAULT_END = "2025-09-14"
SEARCH_ORG = "H5N1"  # search organism/subtype. You may change (H5Nx) if desired.
CLADE_TAG = "2.3.4.4b"  # user-target clade
ENTREZ_TOOL = "agent_h5_clade_agent"
MAFFT_CMD = "mafft"  # must be in PATH

# Cleavage site extraction parameters
# We'll search for the canonical HA cleavage region by looking for the conserved "GLF" motif
# near the HA1/HA2 boundary and extract a window around it. The window size can be changed.
AA_WINDOW = 30  # amino acids before and after the cleavage site to extract
MIN_SEQ_LEN = 1000  # minimal nt length for a full-length HA (heuristic)

# --------------------------- Helper functions ------------------------------

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def fetch_genbank_ha(email, start_date, end_date, out_fasta):
    """Fetch candidate H5 HA sequences from GenBank using Entrez nucleotide search.

    This function performs an Entrez query for HA segment sequences of H5 viruses
    collected within the date range. Beware: GenBank metadata completeness varies.
    """
    Entrez.email = email
    query = f"(" + SEARCH_ORG + f"[Definition] OR H5[All Fields]) AND HA[All Fields] AND \"{start_date}\"[Collection Date] : \"{end_date}\"[Collection Date]"
    print(f"Entrez query: {query}")

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get('IdList', [])
    print(f"Found {len(ids)} candidate GenBank records (ids).\nFetching sequences...")

    if not ids:
        print("No sequences found in GenBank with this query. Consider using GISAID or relaxing the query.")
        return 0

    fetched = 0
    with Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text") as ef:
        seqs = list(SeqIO.parse(ef, "fasta"))
    # Save raw sequences
    SeqIO.write(seqs, out_fasta, "fasta")
    fetched = len(seqs)
    print(f"Wrote {fetched} sequences to {out_fasta}")
    return fetched


def load_fasta_seqs(path):
    if not os.path.exists(path):
        return []
    return list(SeqIO.parse(path, "fasta"))


def deduplicate_seqs(records):
    by_seq = {}
    for r in records:
        seqstr = str(r.seq).upper()
        if seqstr not in by_seq:
            by_seq[seqstr] = r
    return list(by_seq.values())


from Bio.Seq import Seq
def find_cleavage_site_region(seq):
    """
    Find the cleavage site region in an HA nucleotide sequence.
    Returns: (nt_start, nt_end, aa_window, rel_idx) or None
    """
    # 1. Clean sequence
    clean_seq = seq.replace("-", "").replace("n", "N").upper()

    # 2. Pad to multiple of 3
    pad_len = (3 - (len(clean_seq) % 3)) % 3
    if pad_len > 0:
        clean_seq = clean_seq + ("N" * pad_len)

    # 3. Try all three frames
    for frame in range(3):
        try:
            aa = str(Seq(clean_seq[frame:]).translate(to_stop=False))
        except Exception as e:
            print(f"⚠️ Translation error in frame {frame}: {e}")
            continue

        # cleavage motifs to search
        motifs = ["RERRRKR", "RRKR", "RKRK"]

        for motif in motifs:
            if motif in aa:
                pos = aa.find(motif)
                start = max(0, pos - 10)
                end = pos + len(motif) + 10
                aa_window = aa[start:end]

                # nt 좌표 계산 (frame 고려)
                nt_start = frame + start * 3
                nt_end = frame + end * 3
                rel_idx = pos

                return (nt_start, nt_end, aa_window, rel_idx)

    # 못 찾으면 None
    return None

def extract_cleavage_regions_from_fasta(in_fasta, out_fasta, summary_csv):
    records = load_fasta_seqs(in_fasta)
    out_records = []
    rows = []
    for r in records:
        seqstr = str(r.seq)
        if len(seqstr) < MIN_SEQ_LEN:
            # skip short partials
            continue
        found = find_cleavage_site_region(seqstr)
        if found is None:
            rows.append({'accession': r.id, 'found': False, 'note': 'GLF not found'})
            continue
        nt_start, nt_end, aa_window, rel_idx = found
        sub_nt = seqstr[nt_start:nt_end]
        # create a record with header containing origin info if available
        header = r.id
        rec = SeqRecord(Seq(sub_nt), id=header, description=f"cleavage_window:{nt_start}-{nt_end}")
        out_records.append(rec)
        # determine motif (the 1-7 aa immediately preceding GLF)
        # locate GLF in AA window
        aa_full = str(Seq(sub_nt).translate())
        pos = aa_full.find('GLF')
        motif = ''
        if pos != -1:
            # take up to 10 aa before GLF
            motif = aa_full[max(0, pos-10):pos+3]
        rows.append({'accession': r.id, 'found': True, 'motif': motif, 'nt_start': nt_start, 'nt_end': nt_end})
    SeqIO.write(out_records, out_fasta, "fasta")
    df = pd.DataFrame(rows)
    df.to_csv(summary_csv, index=False)
    print(f"Extracted {len(out_records)} cleavage-region sequences to {out_fasta} and wrote summary to {summary_csv}")
    return len(out_records)


def run_mafft(in_fasta, out_fasta, is_protein=False):
    cmd = [MAFFT_CMD, '--auto', in_fasta]
    print('Running MAFFT: ' + ' '.join(cmd))
    with open(out_fasta, 'w') as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        print('MAFFT failed:', p.stderr)
        return False
    return True


def summarize_cleavage_motifs(summary_csv, out_table):
    df = pd.read_csv(summary_csv)
    # collapse motif counts
    motif_counts = df[df['found']==True].groupby('motif').size().reset_index(name='count')
    motif_counts.to_csv(out_table, index=False)
    print(f"Wrote motif summary to {out_table}")
    return motif_counts


# --------------------------- Main agent flow --------------------------------

def main():
    parser = argparse.ArgumentParser(description='Agent pipeline for HA cleavage-site alignment')
    parser.add_argument('--email', required=True, help='Email for Entrez (NCBI)')
    parser.add_argument('--start', default=DEFAULT_START)
    parser.add_argument('--end', default=DEFAULT_END)
    parser.add_argument('--outdir', default='results')
    parser.add_argument('--gisaid_fasta', default=None, help='Optional pre-downloaded GISAID FASTA of HA sequences')
    args = parser.parse_args()

    outdir = args.outdir
    ensure_dir(outdir)
    data_dir = os.path.join(outdir, 'data')
    ensure_dir(data_dir)

    genbank_fasta = os.path.join(data_dir, 'raw_genbank_ha.fasta')
    fetched = fetch_genbank_ha(args.email, args.start, args.end, genbank_fasta)

    combined = []
    combined.extend(load_fasta_seqs(genbank_fasta))
    if args.gisaid_fasta and os.path.exists(args.gisaid_fasta):
        combined.extend(load_fasta_seqs(args.gisaid_fasta))
    combined = deduplicate_seqs(combined)
    combined_fasta = os.path.join(data_dir, 'combined_ha.fasta')
    SeqIO.write(combined, combined_fasta, 'fasta')
    print(f"Wrote {len(combined)} unique combined sequences to {combined_fasta}")

    cleavage_fasta = os.path.join(data_dir, 'ha_cleavage_regions.fasta')
    summary_csv = os.path.join(outdir, 'cleavage_summary.csv')
    n = extract_cleavage_regions_from_fasta(combined_fasta, cleavage_fasta, summary_csv)

    # translate cleavage fasta to protein fasta for alignment
    prot_fasta = os.path.join(data_dir, 'ha_cleavage_regions_prot.fasta')
    prot_recs = []
    for r in load_fasta_seqs(cleavage_fasta):
        aa = str(Seq(str(r.seq)).translate())
        prot_recs.append(SeqRecord(Seq(aa), id=r.id, description=''))
    SeqIO.write(prot_recs, prot_fasta, 'fasta')

    aligned_out = os.path.join(outdir, 'ha_cleavage_alignment.fasta')
    ok = run_mafft(prot_fasta, aligned_out)
    if not ok:
        print('Alignment failed; please ensure MAFFT is installed and accessible in PATH.')

    motif_table = os.path.join(outdir, 'cleavage_motif_table.csv')
    summarize_cleavage_motifs(summary_csv, motif_table)

    print('\nPipeline finished. Key outputs:')
    print(' - Combined sequences: ' + combined_fasta)
    print(' - Cleavage-region sequences: ' + cleavage_fasta)
    print(' - Amino-acid alignment: ' + aligned_out)
    print(' - Motif table: ' + motif_table)
    print('\nCaveats: for complete clade 2.3.4.4b sets, download GISAID EpiFlu data (requires credentials) or use curated Nextstrain builds as a starting point.\n')


if __name__ == '__main__':
    main()
