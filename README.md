# BPG-SH: Building Pan-Genome Based Sequence Homology

## Description

BPG-SH is a bioinformatics pipeline designed to construct pan-genomes and analyze sequence homology. The pipeline involves several steps, from environment setup to gene classification, utilizing various tools and scripts for genome annotation, clustering, and analysis.

## Pipeline Steps

### 1. Activate the Conda Environment

First, you need to activate the Conda environment specified in the `pangenome_env.yml` file. This environment includes all the necessary dependencies for running the pipeline.

```bash
conda env create -f pangenome_env.yml
conda activate pangenome_env
```

### 2. Load Genomes
The 'genome_loader.py' script is used to load genomes. This script will be referenced by other scripts throughout the pipeline.

```bash
python genome_loader.py
```
### 3. Structural Annotation
Step 3.1: Prokka Annotation
Use the 'annotation.py' script to perform genome annotation with Prokka.

```bash
python annotation.py
```
### Step 3.2: Move GFF Files
After running Prokka, manually move the generated GFF files to the gff directory.

```bash
mv path_to_prokka_output/*.gff gff/
```

### Step 3.3: Extract Proteins
Extract proteins from the GFF files using the 'extract_proteins.py' script.

```bash
python extract_proteins.py
```
### Step 3.4: CD-HIT Clustering
Cluster the extracted proteins using CD-HIT through the 'cdhit_clustering.py' script.

```bash
python cdhit_clustering.py
```
### Step 3.5: Parse CD-HIT Results
Parse the results from the CD-HIT clustering using the 'parse_cdhit.py' script.

```bash
python parse_cdhit.py
```
### Step 3.6: Gene Classification
Classify genes into core, dispensable, and unique categories using the 'gene_classifier.py' script.

```bash
python gene_classifier.py
```
