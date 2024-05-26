import os
import subprocess
import logging
import multiprocessing
from genome_loader import load_genomes

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_prokka(fasta_path, output_dir):
    """
    Spustí Prokku k anotaci genomu.

    Parametry:
    fasta_path (str): Cesta ke vstupnímu souboru FASTA.
    output_dir (str): Adresář pro uložení výstupu Prokka.
    """
    base_name = os.path.basename(fasta_path).rsplit('.', 1)[0]
    prokka_output_dir = os.path.join(output_dir, base_name)
    os.makedirs(prokka_output_dir, exist_ok=True)

    # Kontrola počtu dostupných procesorů
    num_cores = multiprocessing.cpu_count()
    requested_cores = min(num_cores, 2)  # Použit 2 cores, nebo méně, pokud je jich k dispozici méně

    prokka_cmd = [
        'prokka', 
        '--outdir', prokka_output_dir, 
        '--prefix', base_name,
        '--genus', 'Escherichia', 
        '--species', 'coli', 
        '--strain', base_name, 
        '--cpus', str(requested_cores), 
        '--force',  # přepsání existujících souborů
        fasta_path
    ]

    result = subprocess.run(prokka_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        logging.info(f"Prokka úspěšně dokončil anotaci pro {fasta_path}")
    else:
        logging.error(f"Prokka se nepodařilo anotovat {fasta_path}: {result.stderr}")

def process_fasta(input_dir, output_dir):
    """
    Zpracuje všechny soubory FASTA ve vstupním adresáři pomocí programu Prokka.

    Parametry:
    input_dir (str): Adresář obsahující vstupní soubory FASTA.
    output_dir (str): Adresář pro uložení výstupů Prokka.
    """
    genomes = load_genomes(input_dir)
    
    for genome_id, record in genomes.items():
        fasta_path = os.path.join(input_dir, f"{genome_id}.fasta")
        if not os.path.isfile(fasta_path):
            logging.error(f"Soubor {fasta_path} neexistuje.")
            continue
        if os.path.getsize(fasta_path) == 0:
            logging.error(f"Soubor {fasta_path} je prázdný.")
            continue
        
        run_prokka(fasta_path, output_dir)

def main():
    """
    Hlavní funkce pro zpracování FASTA souborů a anotaci genomů pomocí Prokka.
    """
    input_dir = 'cesta/k/fasta/souborům' 
    output_dir = 'cesta/do/výstupního/souboru/roary'
    
    process_fasta(input_dir, output_dir)

if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()
