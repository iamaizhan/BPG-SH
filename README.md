# BPG-SH: Building Pan-Genome Based Sequence Homology

## Description

BPG-SH is a bioinformatics pipeline designed to construct pan-genomes and analyze sequence homology. The pipeline involves several steps, from environment setup to gene classification, utilizing various tools and scripts for genome annotation, clustering, and analysis.

**System Architecture**:
- `genome_loader.py`: Loads genome data.
- `annotation.py`: Annotates genome data.
- `extract_proteins.py`: Extracts protein sequences.
- `cdhit_clustering.py`: Performs clustering using CD-HIT.
- `parse_cdhit.py`: Parses CD-HIT results.
- `gene_classifier.py`: Classifies genes.
- `core_phylogenetic_tree.py`: Creates phylogenetic trees.

## Requirements

**Hardware**: A powerful processor and sufficient RAM (recommended at least 16 GB).

**Software**: 
- **Operating System**: Linux/Windows/MacOS
- **Dependencies**: Defined in `pangenome_env.yml`.

## Installation and Setup

1. Install Anaconda or Miniconda.
2. Create and activate the environment:
   ```bash
   conda env create -f pangenome_env.yml
   conda activate pangenome_env

## Pipeline Steps

The following steps outline the pipeline for assembling pan-genomes based on sequence homology:

### Step 1: Load Genome Data
The `genome_loader.py` script is used to load genomes in FASTA formats. This script will be referenced by other scripts throughout the pipeline.

```bash
python genome_loader.py
```
### Step 2: Structural Annotation
Use the `annotation.py` script to perform genome annotation with Prokka.

```bash
python annotation.py
```
### Step 3: Move GFF Files
After running Prokka, manually move the generated GFF files to the gff directory.

```bash
mv path_to_prokka_output/*.gff gff/
```

### Step 4: Extract Proteins
Extract proteins from the GFF files using the `extract_proteins.py` script.

```bash
python extract_proteins.py
```
### Step 5: CD-HIT Clustering
Cluster the extracted proteins using CD-HIT through the `cdhit_clustering.py` script. Clustering helps group similar sequences together.

```bash
python cdhit_clustering.py
```
### Step 6: Parse CD-HIT Results
Parse the results from the CD-HIT clustering using the `parse_cdhit.py` script to make the data usable for subsequent steps.

```bash
python parse_cdhit.py
```
### Step 7: Gene Classification
Classify genes into core, dispensable, and unique categories using the `gene_classifier.py` script.

```bash
python gene_classifier.py
```
### Step 8: Gene Classification
Create Phylogenetic Trees
Generate phylogenetic trees to analyze the evolutionary relationships between the genes using the `core_phylogenetic_tree.py` script.
```bash
python core_phylogenetic_tree.py
```

