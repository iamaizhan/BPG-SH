import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from genome_loader import load_genomes
from gff_loader import load_all_gff_files
from gff_processor import extract_and_translate 

def run_cd_hit(input_fasta, output_dir, identity=0.95, coverage=0.9):
    output_fasta = os.path.join(output_dir, os.path.basename(input_fasta).replace('_proteins', '_cdhit_clusters'))
    cmd = [
        'cd-hit',
        '-i', input_fasta,
        '-o', output_fasta,
        '-c', str(identity),
        '-aL', str(coverage),  # alignment coverage for the longer sequence
        '-aS', str(coverage),  # alignment coverage for the shorter sequence
        '-d', '0'  # description in .clstr file, 0 -> full
    ]
    subprocess.run(cmd, check=True)
    print(f"CD-HIT klastrování dokončeno: {output_fasta}")
