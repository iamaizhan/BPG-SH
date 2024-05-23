import os
import logging
import subprocess

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_cd_hit(input_fasta, output_dir, identity=0.95, coverage=0.99):
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
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith('.fasta')]

    for fasta_file in fasta_files:
        run_cd_hit(fasta_file, output_dir, identity, coverage)

if __name__ == "__main__":
    input_fasta_dir = '/cesta_k/proteins'
    output_dir = '/cesta_k_vystupni_slozce/clstr'
    cluster_all_fastas(input_fasta_dir, output_dir)
