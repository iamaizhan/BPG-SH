import os
import pandas as pd
from collections import defaultdict
import csv

def parse_cd_hit_clstr(clstr_file_path):
    clusters = {}
    cluster_id = None
    with open(clstr_file_path, 'r') as file:
        for line in file:
            if line.startswith('>Cluster'):
                cluster_id = line.strip().split()[1]
                clusters[cluster_id] = []
            else:
                gene_id = line.strip().split('>')[1].split('...')[0]
                clusters[cluster_id].append(gene_id)
    return clusters

def parse_multiple_cd_hit_clstr(clstr_directory):
    all_clusters = {}
    for filename in os.listdir(clstr_directory):
        if filename.endswith('.clstr'):
            genome_id = filename.split('_proteins_cdhit_clusters.fasta.clstr')[0]
            clstr_file_path = os.path.join(clstr_directory, filename)
            all_clusters[genome_id] = parse_cd_hit_clstr(clstr_file_path)
    return all_clusters

def create_presence_absence_matrix(parsed_clusters):
    all_gene_names = set()
    for genome, clusters in parsed_clusters.items():
        for genes in clusters.values():
            all_gene_names.update(genes)
    
    genomes = list(parsed_clusters.keys())
    
    presence_absence_matrix = {gene: [0] * len(genomes) for gene in all_gene_names}
    
    for genome_index, genome in enumerate(genomes):
        clusters = parsed_clusters[genome]
        for genes in clusters.values():
            for gene in genes:
                presence_absence_matrix[gene][genome_index] = 1

    presence_absence_matrix_df = pd.DataFrame.from_dict(presence_absence_matrix, orient='index', columns=genomes)
    presence_absence_matrix_df.reset_index(inplace=True)
    presence_absence_matrix_df.rename(columns={'index': 'Gene'}, inplace=True)
    return presence_absence_matrix_df

def main():
    clstr_directory = '/cesta_k/clstr'
    output_file = '/cesta_k_vystupni_slozce/gene_presence_absence.csv'

    parsed_clusters = parse_multiple_cd_hit_clstr(clstr_directory)
    print(f"Parsované {len(parsed_clusters)} genomů ze souborů klastrů CD-HIT.")

    presence_absence_matrix_df = create_presence_absence_matrix(parsed_clusters)
    
    presence_absence_matrix_df.to_csv(output_file, index=False)
    print("Matrice přítomnosti/nepřítomnosti genů byla úspěšně uložena.")

if __name__ == "__main__":
    main()