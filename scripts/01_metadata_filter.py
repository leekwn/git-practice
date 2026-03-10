import pandas as pd
#Step01
meta = pd.read_excel("gisaid_epiflu_isolates_9990.xls")
print(meta.columns)

#Step02:Standardize text
meta['Location'] = meta['Location'].astype(str)
# Keep only Japan or Korea
selected = meta[meta['Location'].str.contains("Japan|Korea", case=False, na=False)]
print(len(selected))
selected[['Isolate_Id', 'Location']].head()

#Step03
from Bio import SeqIO
# Extract identifiers (simplify to what appears in FASTA header)
selected_ids = set(selected['Isolate_Id'].astype(str).tolist())

input_fasta = "gisaid_epiflu_sequence_8588_isolates.fasta"
output_fasta = "HPAI_JapanKorea_2020_2025.fasta"

count = 0
with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(input_fasta, "fasta"):
        # check if any selected id is in header
        if any(id_ in record.description for id_ in selected_ids):
            SeqIO.write(record, out_f, "fasta")
            count += 1
print(f"{count} sequences written to {output_fasta}")
#
selected.to_csv("metadata_JapanKorea.csv", index=False)


'''
for record in SeqIO.parse("GISAID_HPAI_2020_2025.fasta", "fasta"):
    if "PB2" in record.description:
        segments["PB2"].append(record)
    elif "PB1" in record.description:
        segments["PB1"].append(record)
    elif "PA" in record.description:
        segments["PA"].append(record)
    elif "HA" in record.description:
        segments["HA"].append(record)
    elif "NP" in record.description:
        segments["NP"].append(record)
    elif "NA" in record.description:
        segments["NA"].append(record)
    elif "M" in record.description:
        segments["M"].append(record)
    elif "NS" in record.description:
        segments["NS"].append(record)

for gene, recs in segments.items():
    SeqIO.write(recs, f"{gene}.fasta", "fasta")
'''
