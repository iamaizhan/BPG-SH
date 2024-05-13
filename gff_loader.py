import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genome_loader import load_genomes

def load_gff(gff_filename):
    gff_data = []
    try:
        with open(gff_filename, 'r') as gff_file:
            reader = csv.reader(gff_file, delimiter='\t')
            for row in reader:
                if row[0].startswith('#') or len(row) < 9:
                    continue
                gff_data.append(row)
        return gff_data
    except FileNotFoundError:
        print(f"GFF soubor {gff_filename} nebyl nalezen.")
        return []
    except Exception as e:
        print(f"Chyba při čtení souboru GFF: {e}")
        return []

def load_all_gff_files(directory_path):
    gff_files_data = {}
    try:
        for filename in os.listdir(directory_path):
            if filename.endswith('.gff') or filename.endswith('.gff3'):
                full_path = os.path.join(directory_path, filename)
                gff_data = load_gff(full_path)
                if gff_data:
                    gff_files_data[filename.split('.')[0]] = gff_data 
                    print(f"Načten {filename} s {len(gff_data)} záznamy")
                else:
                    print(f"Žádná data nejsou načtena z {filename}")
    except FileNotFoundError:
        print(f"Adresář {directory_path} nebyl nalezen.")
    except Exception as e:
        print(f"Došlo k neočekávané chybě: {e}")

    return gff_files_data

# directory_path = '/Users/aijan/Desktop/BP/program/gffs'
# all_gff_data = load_all_gff_files(directory_path)
