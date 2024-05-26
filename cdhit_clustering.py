import os
import logging
import subprocess

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_cd_hit(input_fasta, output_dir, identity=0.95, coverage=0.99):
    """
    Spustí CD-HIT ke klastrování proteinových sekvencí.

    Parametry:
    input_fasta (str): Cesta ke vstupnímu FASTA souboru obsahujícímu proteinové sekvence.
    output_dir (str): Adresář pro uložení výstupu CD-HIT.
    identity (float): Práh identity sekvence pro klastrování.
    coverage (float): Práh pokrytí sekvence pro klastrování.
    """
    output_fasta = os.path.join(output_dir, os.path.basename(input_fasta).replace('_proteins', '_cdhit_clusters'))
    cmd = [
        'cd-hit',
        '-i', input_fasta,
        '-o', output_fasta,
        '-c', str(identity),
        '-aL', str(coverage),
        '-aS', str(coverage),
        '-d', '0' 
    ]
    logging.info(f"Spuštění CD-HIT pro {fasta_id} s identitou={identity} a pokrytím={coverage}")
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(f"CD-HIT klastrování bylo úspěšně dokončeno pro: {fasta_id}")
        logging.debug(f"CD-HIT výstup pro {fasta_id}: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"CD-HIT klastrování se nezdařilo pro {fasta_id}: {e.stderr}")

def cluster_all_fastas(fasta_dir, output_dir, identity=0.95, coverage=0.99):
    """
    Shlukuje všechny FASTA soubory v zadaném adresáři pomocí CD-HIT.

    Parametry:
    fasta_dir (str): Cesta ke vstupnímu FASTA souboru obsahujícímu proteinové sekvence.
    output_dir (str): Adresář pro uložení výstupu CD-HIT.
    identity (float): Práh identity sekvence pro klastrování.
    coverage (float): Práh pokrytí sekvence pro klastrování.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fasta_files = []
    for f in os.listdir(fasta_dir):
        if f.endswith('.fasta'):
            fasta_files.append(os.path.join(fasta_dir, f))

    for fasta_file in fasta_files:
        run_cd_hit(fasta_file, output_dir, identity, coverage)

def main():
    """
    Hlavní funkce pro analýzu argumentů a spuštění shlukování CD-HIT.
    """
    input_fasta_dir = '/cesta/k/souborům/proteins'
    output_dir = '/cesta/do/výstupního/souboru/clstr'
    
    cluster_all_fastas(input_fasta_dir, output_dir)
    
if __name__ == "__main__":
    main()
