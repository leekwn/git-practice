# git-practice
git과 github 실습을 위한 저장소
작성일자: 2022.1.15
it changed on 2024.1.06

# Avian Influenza Sequence Analysis Pipeline

A bioinformatics pipeline for processing and analyzing influenza A virus genome sequences, including sequence filtering, segment extraction, alignment, genotype clustering, and comparative seasonal analysis.

This project was developed for the analysis of avian influenza virus datasets (e.g., H5N1, H6N1) obtained from sequence databases such as **GISAID**.

---

# Project Overview

Influenza A viruses contain **eight genomic RNA segments**:

PB2, PB1, PA, HA, NP, NA, M, and NS.

This pipeline performs the following tasks:

1. Metadata filtering and sequence selection
2. Segment extraction from whole genome FASTA files
3. Sequence quality filtering
4. Multiple sequence alignment (MAFFT)
5. Genotype clustering across segments
6. Amino-acid concatenation and variant detection
7. Seasonal comparative analysis (e.g., Korea vs Japan)
8. Visualization and summary statistics

The pipeline helps identify **genotype diversity, evolutionary relationships, and seasonal transmission patterns**.

---

# Repository Structure

```
avian-influenza-sequence-pipeline/

README.md
requirements.txt

scripts/
    01_metadata_filter.py
    02_segment_split.py
    03_sequence_quality_filter.py
    05_genotype_clustering.py
    06_concat_translate_variants.py
    07_season_pairwise_identity_heatmap.py
    08_append_metadata_to_fasta.py

agents/
    h5_clade_cleavage_pipeline.py
    h5n1_literature_agent.py

data/
    raw/
    filtered/
    aligned/

metadata/
```

---

# Pipeline Workflow

The recommended execution order is:

```
1. 01_metadata_filter.py
2. 02_segment_split.py
3. 03_sequence_quality_filter.py
4. MAFFT alignment
5. 05_genotype_clustering.py
6. 06_concat_translate_variants.py
7. 07_season_pairwise_identity_heatmap.py
8. 08_append_metadata_to_fasta.py
```

---

# Step Descriptions

## 1 Metadata Filtering

Extract sequences of interest from metadata and FASTA files (e.g., selecting Japan and Korea isolates).

Input:

* GISAID metadata file
* Raw FASTA sequences

Output:

* Filtered FASTA dataset
* Metadata table

---

## 2 Segment Extraction

Influenza genomes contain eight segments.
This step separates sequences into individual segment FASTA files.

Example output:

```
PB2.fasta
PB1.fasta
PA.fasta
HA.fasta
NP.fasta
NA.fasta
MP.fasta
NS.fasta
```

---

## 3 Sequence Quality Filtering

Removes sequences that:

* are shorter than expected segment length
* contain excessive ambiguous bases
* appear incomplete

This improves alignment quality and downstream analysis.

---

## 4 Multiple Sequence Alignment

Sequences are aligned using **MAFFT**.

Example command:

```
mafft --auto input.fasta > aligned.fasta
```

---

## 5 Genotype Clustering

Aligned sequences are clustered using sequence similarity to define **segment lineages**.

Each isolate receives a genotype based on its segment lineage combination.

Example genotype:

```
PB2-A / PB1-B / PA-A / HA-C / NP-A / NA-B / M-A / NS-A
```

---

## 6 Amino Acid Concatenation and Variant Detection

Steps:

1. Translate nucleotide sequences into amino acids
2. Concatenate proteins across segments
3. Identify variable amino-acid positions

Outputs:

```
concatenated_aa.fasta
aa_variable_sites.csv
```

---

## 7 Seasonal Pairwise Identity Analysis

This analysis compares sequences between **countries and influenza seasons**.

Example comparison:

```
Korea vs Japan
seasonal similarity (HA gene)
```

Outputs:

```
counts_matrix.csv
pct_matrix.csv
heatmap_counts.png
heatmap_pct.png
```

---

## 8 Metadata Annotation

FASTA headers are updated to include metadata fields such as:

```
isolate | country | collection_date | season
```

This allows easier downstream analysis.

---

# Requirements

Python ≥ 3.9

Python packages:

```
biopython
pandas
numpy
scikit-learn
matplotlib
seaborn
tqdm
```

Install dependencies:

```
pip install -r requirements.txt
```

External tools:

* **MAFFT** (for sequence alignment)

Install MAFFT:

```
conda install -c bioconda mafft
```

or

```
sudo apt install mafft
```

---

# Data Sources

Sequence data may be obtained from:

* **GISAID EpiFlu Database**
* **NCBI GenBank**
* other influenza surveillance databases

Please ensure compliance with the **data usage policies** of the respective databases.

---

# Example Use Case

This pipeline can be used to study:

* Influenza genotype evolution
* Cross-regional viral transmission
* Seasonal virus diversity
* Mutation patterns in viral proteins
* HA cleavage site motifs

---

# Reproducibility

To reproduce the analysis:

```
git clone https://github.com/yourusername/avian-influenza-sequence-pipeline.git
cd avian-influenza-sequence-pipeline
pip install -r requirements.txt
```

Place sequence data in:

```
data/raw/
```

Run the pipeline scripts sequentially.

---

# License

This project is provided for research and educational purposes.

---

# Author

Research pipeline developed for influenza virus genomic analysis.

---

# Citation

If you use this pipeline in your research, please cite:

```
Author, Year
Avian Influenza Sequence Analysis Pipeline
GitHub repository
```
