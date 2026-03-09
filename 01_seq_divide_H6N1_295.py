from Bio import SeqIO
from collections import defaultdict

input_fasta = "gisaid_epiflu_sequence_nt_295.fasta"

segments = defaultdict(list)
genes = ["PB2|", "PB1|", "PA|", "HA|", "NP|", "NA|", "MP|", "NS|"]

unmatched = []

for record in SeqIO.parse(input_fasta, "fasta"):
    desc = record.description.upper()
    found = False
    for gene in genes:
        if gene in desc:
            # remove the pipe for clean filenames
            seg = gene.replace("|", "")
            segments[seg].append(record)
            found = True
            break
    if not found:
        unmatched.append(record)

# Write out per-segment FASTA files
for seg, recs in segments.items():
    outfile = f"{seg}_H6N1_world_2010_2025.fasta"
    SeqIO.write(recs, outfile, "fasta")
    print(f"{seg}: {len(recs)} sequences written to {outfile}")

# Optionally write unmatched sequences
if unmatched:
    SeqIO.write(unmatched, "unknown_segment.fasta", "fasta")
    print(f"Unmatched sequences: {len(unmatched)} → written to unknown_segment.fasta")
