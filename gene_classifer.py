import pandas as pd

def classify_genes(presence_absence_matrix):
    total_genomes = presence_absence_matrix.shape[1]
    core_genes = []
    accessory_genes = []
    unique_genes = []

    for gene, row in presence_absence_matrix.iterrows():
        presence_count = row.sum()
        if presence_count == total_genomes:
            core_genes.append(gene)
        elif presence_count == 1:
            unique_genes.append(gene)
        else:
            accessory_genes.append(gene)

    return core_genes, accessory_genes, unique_genes

def save_genes_to_file(genes, filepath):
    try:
        with open(filepath, 'w') as f:
            f.write("\n".join(genes))
        print(f"Geny úspěšně uložené do {filepath}")
    except Exception as e:
        print(f"Chyba při ukládání genů do {filepath}: {e}")

def main():
    input_file = '/cesta_k/gene_presence_absence.csv'
    try:
        presence_absence_matrix = pd.read_csv(input_file, index_col=0)
    except Exception as e:
        print(f"Chyba při čtení vstupního souboru {input_file}: {e}")
        return
    core_genes, accessory_genes, unique_genes = classify_genes(presence_absence_matrix)

    print(f"Core Geny: {len(core_genes)}")
    print(f"Postradatelné Geny: {len(accessory_genes)}")
    print(f"Jedinečné Geny: {len(unique_genes)}")
    
    save_genes_to_file(core_genes, '/cesta_k_vystupni_slozce/core_genes.txt')
    save_genes_to_file(accessory_genes, '/cesta_k_vystupni_slozce/accessory_genes.txt')
    save_genes_to_file(unique_genes, '/cesta_k_vystupni_slozce/unique_genes.txt')

if __name__ == "__main__":
    main()
