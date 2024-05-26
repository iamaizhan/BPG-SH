import os
import pandas as pd
import csv
import seaborn as sns
import matplotlib.pyplot as plt

def parse_cd_hit_clstr(clstr_file_path):
    """
    Analyzuje soubor CD-HIT .clstr a extrahuje klastry.

    Parametry:
    clstr_file_path (str): Cesta k souboru CD-HIT .clstr.

    Vrátí:
    dict: Parsované clustery.
    """
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
    """
    Zpracuje více souborů CD-HIT .clstr z adresáře.

    Parametry:
    clstr_directory (str): Adresář obsahující soubory CD-HIT .clstr.

    Vrátí:
    dict: Parsované klastry pro všechny genomy.
    """
    all_clusters = {}
    for filename in os.listdir(clstr_directory):
        if filename.endswith('.clstr'):
            genome_id = filename.split('_proteins_cdhit_clusters.fasta.clstr')[0]
            clstr_file_path = os.path.join(clstr_directory, filename)
            all_clusters[genome_id] = parse_cd_hit_clstr(clstr_file_path)
    return all_clusters

def create_presence_absence_matrix(parsed_clusters):
    """
    Vytvoří matici přítomnosti/absence z analyzovaných shluků CD-HIT.

    Parametry:
    parsed_clusters (dict): Parsované shluky ze souborů CD-HIT .clstr.

    Vrátí:
    DataFrame: Matice přítomnosti/absence.
    """
    all_gene_names = set()
    for genome, clusters in parsed_clusters.items():
        for genes in clusters.values():
            all_gene_names.update(genes)
    
    genomes = list(parsed_clusters.keys())
    
    presence_absence_matrix = {}
    for gene in all_gene_names:
        presence_absence_matrix[gene] = [0] * len(genomes)
    
    for genome_index, genome in enumerate(genomes):
        clusters = parsed_clusters[genome]
        for genes in clusters.values():
            for gene in genes:
                presence_absence_matrix[gene][genome_index] = 1

    presence_absence_matrix_df = pd.DataFrame.from_dict(presence_absence_matrix, orient='index', columns=genomes)
    presence_absence_matrix_df.reset_index(inplace=True)
    presence_absence_matrix_df.rename(columns={'index': 'Gene'}, inplace=True)
    return presence_absence_matrix_df

def heatmap(input_file, output_image):
    """
    Vykreslí setříděnou heatmapu matice přítomnosti/absence

    Parametry:
    input_file (str): Cesta k souboru CSV obsahujícímu matici přítomnosti/absence.
    output_image (str): Cesta k uložení výstupního obrázku heatmapy.
    """
    presence_absence_matrix = pd.read_csv(input_file, index_col=0)
    
    presence_absence_matrix = presence_absence_matrix.loc[presence_absence_matrix.sum(axis=1).sort_values(ascending=False).index]
    presence_absence_matrix = presence_absence_matrix[presence_absence_matrix.sum().sort_values(ascending=True).index]

    plt.figure(figsize=(15, 10))
    sns.heatmap(presence_absence_matrix, cmap="viridis", cbar=False)
    plt.title('Mapa přítomnosti/absence genů')
    plt.xlabel('Genomy')
    plt.ylabel('Geny')
    plt.savefig(output_image)
    plt.show()

def main():
    """
    Hlavní funkce pro zpracování shluků CD-HIT a vykreslení tepelné mapy přítomnosti/absence.
    """
    clstr_directory = '/cesta/k/souborům/clstr'
    output_file = '/cesta/do/output/gene_presence_absence.csv'

    parsed_clusters = parse_multiple_cd_hit_clstr(clstr_directory)
    print(f"Parsované {len(parsed_clusters)} genomů ze souborů klastrů CD-HIT.")

    presence_absence_matrix_df = create_presence_absence_matrix(parsed_clusters)
    
    presence_absence_matrix_df.to_csv(output_file, index=False)
    print("Matrice přítomnosti/nepřítomnosti genů byla úspěšně uložena.")

    heatmap(output_file, heatmap_output)
    
if __name__ == "__main__":
    main()
