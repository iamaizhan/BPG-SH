import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def load_genomes(directory_path):
    genome_dict = {}
    try:
        all_files = os.listdir(directory_path)
        fasta_files = [file for file in all_files if file.endswith('.fasta')]
        for file in fasta_files:
            path = os.path.join(directory_path, file)
            genome_id = file.split('.')[0]
            genome_dict[genome_id] = SeqIO.read(path, "fasta")
        print(f"Soubory byly úspěšně načteny: {len(fasta_files)} souborů.")
        for file in fasta_files:
            print(file)
    except FileNotFoundError:
        print(f"Adresář {directory_path} nebyl nalezen.")
    except Exception as e:
        print(f"Došlo k neočekávané chybě: {e}")
    return genome_dict

if __name__ == "__main__":
    path_to_genomes = '/cesta/k/e.coli' # Cesta do .fasta souborů
    load_genomes(path_to_genomes)
