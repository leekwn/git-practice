#!/usr/bin/env python3
"""
concat_translate_variants.py

Usage:
  python concat_translate_variants.py --genes PB2 PB1 PA HA NP NA M NS \
    --indir ./genes_fasta --outdir ./results --msa_program mafft

Outputs:
  - results/<gene>_aligned_nt.fasta   (MAFFT nucleotide alignment per gene)
  - results/<gene>_aligned_aa.fasta   (translated AA alignment per gene)
  - results/concatenated_aa.fasta     (concatenated AA sequences)
  - results/aa_variable_sites.csv     (table of variable AA positions)
"""

import os
import sys
import argparse
import subprocess
from collections import defaultdict, OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ----------------------------- helpers -----------------------------
def run_mafft(input_fasta, output_fasta):
    """Run MAFFT --auto; requires mafft in PATH."""
    cmd = ["mafft", "--auto", "--thread", "1", input_fasta]
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        with open(output_fasta, "w") as fh:
            fh.write(p.stdout)
        return True
    except Exception as e:
        print("MAFFT failed or not installed:", e, file=sys.stderr)
        return False

def read_fasta_to_dict(path):
    seqs = OrderedDict()
    for r in SeqIO.parse(path, "fasta"):
        seqs[r.id] = str(r.seq).upper().replace("U","T")
    return seqs

def write_fasta_from_dict(seqdict, path, description=""):
    recs = [SeqRecord(Seq(s), id=rid, description="") for rid, s in seqdict.items()]
    SeqIO.write(recs, path, "fasta")

def translate_aligned_nt_to_aa(nt_seq):
    """Translate an aligned nucleotide sequence (may include '-') to an AA sequence.
    Gaps that break a codon -> produce '-' for that AA position.
    """
    seq = nt_seq.upper()
    aa_chars = []
    # process in codons
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            # incomplete codon at end -> treat as gap
            aa_chars.append('-')
            continue
        if '-' in codon:
            aa_chars.append('-')
        else:
            aa = str(Seq(codon).translate(table=1))  # standard table
            aa_chars.append(aa)
    return "".join(aa_chars)

# ----------------------------- main pipeline -----------------------------
def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    gene_aa_alignments = OrderedDict()
    gene_len_aa = OrderedDict()
    sample_ids = None

    for gene in args.genes:
        inp = os.path.join(args.indir, f"{gene}.fasta")
        if not os.path.exists(inp):
            print(f"Missing input for gene {gene}: {inp}", file=sys.stderr)
            sys.exit(1)

        # Read raw sequences (unaligned nucleotides)
        seqs = read_fasta_to_dict(inp)
        if sample_ids is None:
            sample_ids = list(seqs.keys())
        else:
            # ensure IDs consistent
            if set(sample_ids) != set(seqs.keys()):
                print(f"Sequence IDs in {gene}.fasta don't match previous genes.", file=sys.stderr)
                print("IDs in first gene:", sample_ids, file=sys.stderr)
                print("IDs in this gene:", list(seqs.keys()), file=sys.stderr)
                sys.exit(1)

        # write a temp fasta for MAFFT
        tmp_in = os.path.join(args.outdir, f"{gene}._tmp_in.fasta")
        write_fasta_from_dict(seqs, tmp_in)

        aligned_nt = os.path.join(args.outdir, f"{gene}_aligned_nt.fasta")
        success = False
        if args.msa_program and args.msa_program.lower() == "mafft":
            success = run_mafft(tmp_in, aligned_nt)
        if not success:
            # fallback: write unaligned (NOT IDEAL)
            print(f"Falling back to unaligned sequences for gene {gene} (NOT recommended).", file=sys.stderr)
            write_fasta_from_dict(seqs, aligned_nt)

        # read aligned nt sequences
        aligned_seqs = read_fasta_to_dict(aligned_nt)

        # translate each aligned nt -> aa alignment
        aa_seqs = OrderedDict()
        for rid in sample_ids:
            aa = translate_aligned_nt_to_aa(aligned_seqs[rid])
            aa_seqs[rid] = aa

        # sanity: all aa seqs must have same length
        lengths = {len(x) for x in aa_seqs.values()}
        if len(lengths) != 1:
            print(f"Inconsistent AA alignment lengths for gene {gene}: {lengths}", file=sys.stderr)
            sys.exit(1)
        gene_len_aa[gene] = lengths.pop()
        gene_aa_alignments[gene] = aa_seqs

        # write AA alignment per gene
        out_aa = os.path.join(args.outdir, f"{gene}_aligned_aa.fasta")
        write_fasta_from_dict(aa_seqs, out_aa)
        # cleanup tmp
        try:
            os.remove(tmp_in)
        except:
            pass

    # Concatenate genes in order (AA)
    concat_aa = OrderedDict((rid, "") for rid in sample_ids)
    gene_boundaries = []  # list of (gene, start_pos_1based, end_pos_1based)
    pos_cursor = 1
    for gene in args.genes:
        L = gene_len_aa[gene]
        for rid in sample_ids:
            concat_aa[rid] += gene_aa_alignments[gene][rid]
        gene_boundaries.append((gene, pos_cursor, pos_cursor + L - 1))
        pos_cursor += L

    # write concatenated AA
    concat_path = os.path.join(args.outdir, "concatenated_aa.fasta")
    write_fasta_from_dict(concat_aa, concat_path)

    # find variable sites in concatenated AA alignment
    n_total = len(next(iter(concat_aa.values())))
    var_sites = []
    for i in range(n_total):
        col = [concat_aa[rid][i] for rid in sample_ids]
        unique = sorted(set([c for c in col if c != '-']))
        if len(unique) > 1:
            # variable site
            # map to gene
            pos1 = i + 1
            gene_name = None
            gene_pos = None
            for (g, start, end) in gene_boundaries:
                if start <= pos1 <= end:
                    gene_name = g
                    gene_pos = pos1 - start + 1
                    break
            ref_aa = concat_aa[sample_ids[0]][i]  # use first sequence as reference
            entry = {
                "concat_pos": pos1,
                "gene": gene_name,
                "gene_aa_pos": gene_pos,
                "ref_aa": ref_aa
            }
            for rid in sample_ids:
                entry[rid] = concat_aa[rid][i]
            var_sites.append(entry)

    # write CSV of variable sites
    import csv
    csv_path = os.path.join(args.outdir, "aa_variable_sites.csv")
    with open(csv_path, "w", newline="") as fh:
        fieldnames = ["concat_pos", "gene", "gene_aa_pos", "ref_aa"] + sample_ids
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for e in var_sites:
            writer.writerow(e)

    print("Done.")
    print("Outputs written to:", args.outdir)
    print(" - per-gene aligned nt: <gene>_aligned_nt.fasta")
    print(" - per-gene aligned aa: <gene>_aligned_aa.fasta")
    print(" - concatenated aa fasta:", os.path.basename(concat_path))
    print(" - variable sites csv:", os.path.basename(csv_path))

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Concatenate 8 genes, translate, and extract AA differences.")
    p.add_argument("--genes", nargs="+", required=True, help="List of genes in desired concatenation order")
    p.add_argument("--indir", required=True, help="Input directory with <gene>.fasta files")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--msa_program", default="mafft", help="MSA program (mafft recommended)")
    args = p.parse_args()
    main(args)
