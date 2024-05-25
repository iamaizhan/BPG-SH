import os
import csv
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genome_loader import load_genomes

# Nastavení základního logování pro záznam průběhu zpracování
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extract_cds_from_gff(gff_file, genome_sequence):
    proteins = []
    try:
        with open(gff_file, 'r') as gff_handle:
            reader = csv.reader(gff_handle, delimiter='\t')
            for row in reader:
                if row[0].startswith('#') or len(row) < 9:
                    continue
                seq_id, source, feature_type, start, end, score, strand, phase, attributes = row
                if feature_type == 'CDS':
                    protein = parse_protein_info(start, end, strand, attributes, genome_sequence)
                    if protein:
                        proteins.append(protein)
    except Exception as e:
        logging.error(f"Nepodařilo se zpracovat GFF soubor {gff_file}: {e}")
    return proteins

def parse_protein_info(start, end, strand, attributes, genome_sequence):
    attributes_dict = {attr.split('=')[0]: attr.split('=')[1] for attr in attributes.split(';') if '=' in attr}
    protein_id = attributes_dict.get('ID', 'unknown_protein')
    gene_name = attributes_dict.get('gene', 'unknown_gene')
    product = attributes_dict.get('product', 'unknown_product')
    start = int(start) - 1
    end = int(end)
    if strand == '+':
        nuc_seq = genome_sequence.seq[start:end]
    else:
        nuc_seq = genome_sequence.seq[start:end].reverse_complement()
    
    try:
        protein_seq = nuc_seq.translate(to_stop=True)
        if not protein_seq:
            raise ValueError("Prázdná sekvence proteinu")
        return SeqRecord(protein_seq, id=gene_name, description=product)
    except Exception as e:
        logging.error(f"Chyba při translace sekvence pro {protein_id} (product: {product}): {e}")
        return None

def save_proteins_to_fasta(proteins, output_file):
    valid_proteins = [p for p in proteins if p is not None]
    if not valid_proteins:
        logging.warning(f"Nebyly nalezeny žádné validní proteiny. Přeskočení souboru {output_file}.")
        return
    try:
        SeqIO.write(valid_proteins, output_file, "fasta")
        logging.info(f"Proteiny uložené do {output_file}")
    except Exception as e:
        logging.error(f"Nepodařilo se zapsat proteiny do {output_file}: {e}")

def process_gffs_and_genomes(gff_directory, genome_directory, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    genome_sequences = load_genomes(genome_directory)
    
    for filename in os.listdir(gff_directory):
        if filename.endswith('.gff') or filename.endswith('.gff3'):
            gff_path = os.path.join(gff_directory, filename)
            genome_id = filename.split('.')[0]
            if genome_id in genome_sequences:
                proteins = extract_cds_from_gff(gff_path, genome_sequences[genome_id])
                output_file = os.path.join(output_directory, f"{genome_id}_proteins.fasta")
                save_proteins_to_fasta(proteins, output_file)
            else:
                logging.error(f"Nebyla nalezena žádná sekvence genomu {genome_id}")


if __name__ == "__main__":
    gff_directory = '/cesta_k/roary/gff'
    genome_directory = '/cesta_k_genomu'
    output_directory = '/cesta_k_vystupni_slozce/proteins'
    
    process_gffs_and_genomes(gff_directory, genome_directory, output_directory)
