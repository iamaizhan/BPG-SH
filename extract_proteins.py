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
    """
    Extrahuje proteinové sekvence ze souboru GFF a zapíše je do souboru FASTA.

    Parametry:
    gff_file (str): Cesta ke vstupnímu GFF souboru.
    genome_sequence (str): Cesta k výstupnímu FASTA souboru.

    Vratí:
    list: Seznam SeqRecords obsahující proteinové sekvence.
    """
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
    """
    Vypíše informace o proteinech z atributů GFF a sekvence genomu.

    Parametry:
    start (str): Počáteční pozice CDS.
    end (str): Koncová pozice CDS.
    strand (str): Informace o vlákně ("+" nebo "-").
    attributes (str): Řetězec atributů ze souboru GFF.
    genome_sequence (SeqRecord): Záznam sekvence genomu.

    Vrátí:
    SeqRecord: SeqRecord přeložené sekvence proteinu.
    """
    attributes_dict = {}
    for attr in attributes.split(';'):
        if '=' in attr:
            key, value = attr.split('=')
            attributes_dict[key] = value
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
    """
    Uloží proteinové sekvence do FASTA souboru.

    Parametry:
    proteins (list): Seznam SeqRecords obsahující proteinové sekvence.
    output_file (str): Cesta k výstupnímu souboru 
    """
    valid_proteins = []
    for p in proteins:
        if p is not None:
            valid_proteins.append(p)
    if not valid_proteins:
        logging.warning(f"Nebyly nalezeny žádné validní proteiny. Přeskočení souboru {output_file}.")
        return
    try:
        SeqIO.write(valid_proteins, output_file, "fasta")
        logging.info(f"Proteiny uložené do {output_file}")
    except Exception as e:
        logging.error(f"Nepodařilo se zapsat proteiny do {output_file}: {e}")

def process_gffs_and_genomes(gff_directory, genome_directory, output_directory):
    """
    Zpracuje všechny soubory GFF ve vstupním adresáři za účelem extrakce proteinů a uloží je do souborů FASTA.
    
    Parametry:
    gff_directory (str): Adresář obsahující vstupní GFF soubory.
    genome_directory (str): Adresář obsahující sekvence genomu.
    output_directory (str): Adresář pro uložení výstupních FASTA souborů.
    """
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

def main():
    """
    Hlavní funkce pro zpracování souborů GFF a sekvencí genomu, extrakci proteinů a jejich uložení do FASTA souborů.
    """
    gff_directory = '/cesta/k/souborům/roary/gff'
    genome_directory = '/cesta/k/fasta/souborům'
    output_directory = '/cesta/do/výstupního/souboru/proteins'
    
    process_gffs_and_genomes(gff_directory, genome_directory, output_directory)

if __name__ == "__main__":
    main()
