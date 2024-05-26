import pandas as pd
import matplotlib.pyplot as plt

def classify_genes(presence_absence_matrix):
    """
    Klasifikuje geny na základě matice přítomnosti/absence.

    Parametry:
    presence_absence_matrix_file (str): Matice přítomnosti/absence.

    Vrátí:
    listy: Seznamy klasifikovaných genů na core, postradatelné a jedinečné.
    """
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
    """
    Uloží klasifikované geny do souboru.

    Parametry:
    genes (list): Seznam genů k uložení.
    filepath (str): Cesta k výstupnímu souboru.
    """
    try:
        with open(filepath, 'w') as f:
            f.write("\n".join(genes))
        print(f"Geny úspěšně uložené do {filepath}")
    except Exception as e:
        print(f"Chyba při ukládání genů do {filepath}: {e}")

def main():
    """
    Hlavní funkce pro klasifikaci genů, uložení výsledků a vykreslení rozložení genů.
    """
    input_file = '/cesta/k/output/gene_presence_absence.csv'
    try:
        presence_absence_matrix = pd.read_csv(input_file, index_col=0)
    except Exception as e:
        print(f"Chyba při čtení vstupního souboru {input_file}: {e}")
        return
    core_genes, accessory_genes, unique_genes = classify_genes(presence_absence_matrix)

    print(f"Core Geny: {len(core_genes)}")
    print(f"Postradatelné Geny: {len(accessory_genes)}")
    print(f"Jedinečné Geny: {len(unique_genes)}")
    
    save_genes_to_file(core_genes, '/cesta/do/output/core_genes.txt')
    save_genes_to_file(accessory_genes, '/cesta/do/output/accessory_genes.txt')
    save_genes_to_file(unique_genes, '/cesta/do/output/unique_genes.txt')

    # Koláčový graf
    labels = ['Core Geny', 'Postradatelné Geny', 'Jedinečné Geny']
    sizes = [len(core_genes), len(accessory_genes), len(unique_genes)]
    colors = ['gold', 'lightcoral', 'lightskyblue']
    explode = (0.1, 0, 0)
    
    plt.figure(figsize=(10, 8))
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%',
            shadow=True, startangle=140)
    plt.axis('equal')  
    plt.title('Distribuce Core, Postradatelých a Jedinečných Genů')
    plt.savefig('gene_classification_pie_chart.png')
    plt.show()

if __name__ == "__main__":
    main()
