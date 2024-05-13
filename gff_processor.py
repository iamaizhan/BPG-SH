import csv
import os
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def load_and_process_gff(directory_path):
    all_proteins = {}
    for filename in os.listdir(directory_path):
        if filename.endswith('.gff') or filename.endswith('.gff3'):
            gff_path = os.path.join(directory_path, filename)
            genome_id = filename.split('.')[0]
            proteins = extract_proteins_from_gff(gff_path, genome_id)
            all_proteins[genome_id] = proteins
    return all_proteins

def extract_proteins_from_gff(gff_path, genome_id):
    proteins = []
    try:
        with open(gff_path, 'r') as gff_file:
            reader = csv.reader(gff_file, delimiter='\t')
            for row in reader:
                if row[0].startsWith('#') or len(row) < 9:
                    continue
                seq_id, source, feature_type, start, end, score, strand, phase, attributes = row
                if feature_type == 'CDS':
                    proteins.append(parse_protein_info(seq_id, start, end, strand, phase, attributes))
    except Exception as e:
        logging.error(f"Failed to process GFF file {gff_path}: {e}")
    return proteins

def parse_protein_info(seq_id, start, end, strand, phase, attributes):
    attributes_dict = {attr.split('=')[0]: attr.split('=')[1] for attr in attributes.split(';') if '=' in attr}
    return SeqRecord(Seq(attributes_dict.get('sequence', '')), id=attributes_dict.get('ID', 'unknown_protein'))
