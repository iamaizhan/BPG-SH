import os
from Bio import SeqIO

def load_genomes(directory_path):
    '''
    Načte genomy z FASTA souborů v daném adresáři.

    Parametry:
    directory_path (str): Cesta k adresáři obsahujícímu FASTA soubory.

    Vrátí:
    dict: Slovník, kde klíče jsou ID genomů (odvozené z názvů souborů) a hodnoty jsou sekvence genomů ve formátu Biopython SeqRecord.
    '''
    genome = {}
    try:
        all_files = os.listdir(directory_path)
        fasta_files = [file for file in all_files if file.endswith('.fasta')]
        for file in fasta_files:
            path = os.path.join(directory_path, file)
            genome_id = file.split('.')[0]
            genome[genome_id] = SeqIO.read(path, "fasta")
        print(f"Soubory byly úspěšně načteny: {len(fasta_files)} souborů.")
        for file in fasta_files:
            print(file)
    except FileNotFoundError:
        print(f"Adresář {directory_path} nebyl nalezen.")
    except Exception as e:
        print(f"Došlo k neočekávané chybě: {e}")
    return genome

if __name__ == "__main__":
    path_to_genomes = '/cesta/k/fasta/souborům/'
    load_genomes(path_to_genomes)
