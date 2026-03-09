#!/usr/bin/env python3
"""
Filter out short or ambiguous influenza A sequences (e.g. H6N1) before alignment.

- Applies to all FASTA files in the input folder (one per segment).
- By default uses approximate expected lengths for influenza A segments.
- Optionally auto-estimates expected length per file (median length).

Defaults:
  min_fraction = 0.9   → keep sequences ≥ 90% of expected length
  max_ambig_frac = 0.02 → drop sequences with > 2% ambiguous (non-ATGC) bases
"""

from Bio import SeqIO
import os
import argparse
import statistics

# Approximate nucleotide lengths for influenza A segments
EXPECTED_LENGTHS = {
    "PB2": 2341,
    "PB1": 2341,
    "PA": 2233,
    "HA": 1700,
    "NP": 1565,
    "NA": 1413,
    "M": 1027,   # use M consistently
    "NS": 890
}

def infer_segment_from_filename(basename):
    """
    Infer segment name from filename.
    Expects something like:
      PB2_JapanKorea.fasta  → PB2
      HA_H6N1.fasta         → HA
    """
    seg_guess = basename.split("_")[0].upper()
    # Map MP → M if needed
    if seg_guess == "MP":
        seg_guess = "M"
    return seg_guess

def estimate_expected_length(records):
    """
    Estimate expected length as the median length of sequences in the file.
    Useful for new subtypes (e.g. H6N1) or mixed datasets.
    """
    lengths = [len(str(rec.seq)) for rec in records if len(rec.seq) > 0]
    if not lengths:
        return None
    return int(round(statistics.median(lengths)))

def filter_sequences(input_fasta, output_fasta,
                     min_fraction=0.9,
                     max_ambig_frac=0.02,
                     use_auto_expected=False):
    """Filter out sequences that are too short or have too many ambiguous bases."""
    basename = os.path.basename(input_fasta)
    seg_guess = infer_segment_from_filename(basename)

    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        print(f"{basename}: no records found, skipping.")
        return

    expected_len = EXPECTED_LENGTHS.get(seg_guess, None)

    if use_auto_expected or expected_len is None:
        auto_len = estimate_expected_length(records)
        if auto_len is None:
            print(f"{basename}: could not estimate expected length, skipping.")
            return
        expected_len = auto_len
        print(f"{basename}: using auto-estimated expected length = {expected_len} nt for segment {seg_guess}")
    else:
        print(f"{basename}: using predefined expected length = {expected_len} nt for segment {seg_guess}")

    kept, removed = 0, 0
    kept_records = []

    for rec in records:
        seq = str(rec.seq).upper()
        seq_len = len(seq)

        if seq_len == 0:
            removed += 1
            continue

        # Count ambiguous bases (non-ATGC; includes N, gaps, etc.)
        ambig = sum(base not in "ATGC" for base in seq)
        ambig_frac = ambig / seq_len

        # Length and ambiguity filters
        if seq_len < min_fraction * expected_len:
            removed += 1
            continue
        if ambig_frac > max_ambig_frac:
            removed += 1
            continue

        kept += 1
        kept_records.append(rec)

    if kept_records:
        SeqIO.write(kept_records, output_fasta, "fasta")
    else:
        # Empty file: optionally you can choose not to write it
        open(output_fasta, "w").close()

    print(f"{basename}: kept {kept}, removed {removed}")
    
    

def main():
    parser = argparse.ArgumentParser(description="Filter out short or ambiguous influenza sequences.")
    parser.add_argument("--indir", required=True, help="Directory containing per-segment FASTA files")
    parser.add_argument("--outdir", default="filtered", help="Output directory for filtered FASTAs")
    parser.add_argument("--min_fraction", type=float, default=0.9,
                        help="Minimum fraction of expected length (default: 0.9)")
    parser.add_argument("--max_ambig_frac", type=float, default=0.02,
                        help="Maximum fraction of ambiguous bases (default: 0.02)")
    parser.add_argument("--auto_expected", action="store_true",
                        help="Auto-estimate expected length per file from median sequence length")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for fname in os.listdir(args.indir):
        if not fname.lower().endswith(".fasta"):
            continue
        inpath = os.path.join(args.indir, fname)
        outpath = os.path.join(args.outdir, fname)
        filter_sequences(
            inpath,
            outpath,
            min_fraction=args.min_fraction,
            max_ambig_frac=args.max_ambig_frac,
            use_auto_expected=args.auto_expected
        )

if __name__ == "__main__":
    main()
