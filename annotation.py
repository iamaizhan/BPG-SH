import os
import subprocess
import logging
import multiprocessing
from genome_loader import load_genomes

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_prokka(fasta_path, output_dir):
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

def process_genomes(genomes, output_dir):
    for genome_id, record in genomes.items():
        run_prokka(record, genome_id, output_dir)

if __name__ == "__main__":
    input_dir = '/cesta/k/e.coli' 
    output_dir = '/cesta/k/výstupní/složce/gff'
    
    genomes = load_genomes(input_dir)
    process_genomes(genomes, output_dir)
