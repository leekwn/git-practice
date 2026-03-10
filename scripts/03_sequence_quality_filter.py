#!/usr/bin/env python3
"""
Filter out short or ambiguous sequences before alignment.
Applies to all FASTA files in the input folder (one per segment).
"""

from Bio import SeqIO
import os
import argparse

# Expected lengths for influenza A segments (approximate H5Nx)
EXPECTED_LENGTHS = {
    "PB2": 2341,
    "PB1": 2341,
    "PA": 2233,
    "HA": 1700,
    "NP": 1565,
    "NA": 1413,
    "MP": 1027,
    "M": 1027,
    "NS": 890
}

def filter_sequences(input_fasta, output_fasta, min_fraction=0.95, max_ambig_frac=0.01):
    """Filter out sequences that are too short or have too many ambiguous bases"""
    basename = os.path.basename(input_fasta)
    seg_guess = basename.split("_")[0].upper()  # e.g., PB2_JapanKorea.fasta → PB2

    expected_len = EXPECTED_LENGTHS.get(seg_guess, None)
    if expected_len is None:
        print(f"⚠️  Unknown segment for {basename}, skipping expected length check.")
        expected_len = 1000  # generic fallback

    kept, removed = 0, 0
    kept_records = []

    for rec in SeqIO.parse(input_fasta, "fasta"):
        seq = str(rec.seq).upper()
        seq_len = len(seq)

        # Count ambiguous bases (non-ATGC)
        ambig = sum(base not in "ATGC" for base in seq)
        ambig_frac = ambig / seq_len if seq_len > 0 else 1.0

        # Length and ambiguity filters
        if seq_len < min_fraction * expected_len:
            removed += 1
            continue
        if ambig_frac > max_ambig_frac:
            removed += 1
            continue

        kept += 1
        kept_records.append(rec)

    SeqIO.write(kept_records, output_fasta, "fasta")
    print(f"{basename}: kept {kept}, removed {removed}")

def main():
    parser = argparse.ArgumentParser(description="Filter out short or ambiguous influenza sequences.")
    parser.add_argument("--indir", required=True, help="Directory containing per-segment FASTA files")
    parser.add_argument("--outdir", default="filtered", help="Output directory for filtered FASTAs")
    parser.add_argument("--min_fraction", type=float, default=0.9, help="Minimum fraction of expected length")
    parser.add_argument("--max_ambig_frac", type=float, default=0.02, help="Maximum fraction of ambiguous bases")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for fname in os.listdir(args.indir):
        if not fname.lower().endswith(".fasta"):
            continue
        inpath = os.path.join(args.indir, fname)
        outpath = os.path.join(args.outdir, fname)
        filter_sequences(inpath, outpath, args.min_fraction, args.max_ambig_frac)

if __name__ == "__main__":
    main()
