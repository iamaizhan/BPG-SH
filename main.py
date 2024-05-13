import os 
from genome_loader import load_genomes
from ąnnotation import run_prokka, annotate_all_genomes
from file_manager import gather_gff_files
from gff_processor import load_and_process_gff
from cd_hit_clustering import run_cd_hit
from clustering_analysis import parse_cdhit_clstr_file, categorize_clusters
from Bio import SeqIO

def main():
    input_dir = '/Users/aijan/Desktop/MAIN PIPELINE/e.coli'
    output_dir = '/Users/aijan/Desktop/MAIN PIPELINE/gffs'
    gff_dir = '/Users/aijan/Desktop/MAIN PIPELINE/gffs'

    # Anotace 
    source_dirs = []
    annotate_all_genomes(input_dir, output_dir)
    for filename in os.listdir(output_dir):
        if os.path.isdir(os.path.join(output_dir, filename)):
            source_dirs.append(os.path.join(output_dir, filename))

    # Každý GFF do jedné složky
    gather_gff_files(source_dirs, gff_dir)

    # Načtení a zpracování GFF -> extrakce proteinů
    all_proteins = load_and_process_gff(gff_dir)
    
    # Klastrování pomocí CD-HIT
    for genome_id, proteins in all_proteins.items():
        protein_fasta_path = os.path.join(output_dir, f"{genome_id}_proteins.fasta")
        SeqIO.write(proteins, protein_fasta_path, "fasta")
        run_cd_hit(protein_fasta_path, output_dir)

    clstr_dir = '/Users/aijan/Desktop/MAIN PIPELINE/output'
    total_genomes = 10
    all_clusters = {}

    for file in os.listdir(clstr_dir):
        if file.endswith('.clstr'):
            file_path = os.path.join(clstr_dir, file)
            clusters = parse_cdhit_clstr_file(file_path)
            for key, value in clusters.items():
                if key in all_clusters:
                    all_clusters[key].extend(value)
                else:
                    all_clusters[key] = value

    core, accessory, unique = categorize_clusters(all_clusters, total_genomes)

    print(f"Jádrové Klastry Genomů: {len(core)}")
    print(f"Postradatelné Klastry Genomů: {len(accessory)}")
    print(f"Jedinečné Klastry Genomů: {len(unique)}")

if __name__ == "__main__":
    main()
