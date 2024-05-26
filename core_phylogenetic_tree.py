import os
import subprocess
from Bio import SeqIO
from BCBio import GFF
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

def load_core_genes(core_genes_file):
    """
    Načte core geny ze zadaného souboru.

    Parametry:
    core_genes_file (str): Cesta k souboru obsahujícímu seznam core genů.

    Vrátí:
    list: Seznam core genů.
    """
    with open(core_genes_file, 'r') as file:
        core_genes = []
        for line in file:
            core_genes.append(line.strip())
    print(f"Načteny core geny: {core_genes}")
    return core_genes

def extract_core_sequences(gff_file, fasta_file, core_genes):
    """
    Extrahuje core sekvence ze zadaných souborů GFF a FASTA.

    Parametry:
    gff_file (str): Cesta k souboru GFF.
    fasta_file (str): Cesta k souboru FASTA.
    core_genes (list): Seznam core genů.

    Vrátí:
    dict: Slovník s názvy genů jako klíči a sekvencemi jako hodnotami.
    """
    sequences = {}
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    with open(gff_file) as gff_handle:
        for rec in GFF.parse(gff_handle, base_dict=fasta_sequences):
            for feature in rec.features:
                if feature.type == "CDS" and 'Name' in feature.qualifiers:
                    gene_name = feature.qualifiers['Name'][0]
                    if gene_name in core_genes:
                        seq = feature.extract(rec.seq)
                        sequences[gene_name] = str(seq)
    print(f"Extrahované sekvence pro geny: {list(sequences.keys())} ze souboru {gff_file}")
    return sequences

def concatenate_core_genes(core_genes, gff_directory, fasta_directory, output_file):
    """
    Spojí core geny pro každý genom a uloží je do jednoho FASTA souboru.

    Parametry:
    core_genes (list): Seznam core genů.
    gff_directory (str): Adresář obsahující GFF soubory.
    fasta_directory (str): Adresář obsahující FASTA soubory.
    output_file (str): Cesta k výstupnímu FASTA souboru.
    """
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(gff_directory):
            if filename.endswith('.gff') or filename.endswith('.gff3'):
                genome_id = os.path.splitext(filename)[0]
                concatenated_seq = ''
                gff_path = os.path.join(gff_directory, filename)
                fasta_path = os.path.join(fasta_directory, genome_id + '.fasta')
                print(f"Zpracování genomu {genome_id}")
                sequences = extract_core_sequences(gff_path, fasta_path, core_genes)
                if sequences:
                    print(f"Nalezena sekvence pro genom {genome_id}")
                for gene in core_genes:
                    if gene in sequences:
                        concatenated_seq += sequences[gene]
                    else:
                        print(f"Gen {gene} nebyl nalezen v genomu {genome_id}")
                if concatenated_seq:
                    outfile.write(f'>{genome_id}\n{concatenated_seq}\n')
                else:
                    print(f"Nebyly nalezeny žádné sekvence pro genom {genome_id}")

def run_mafft(input_fasta, output_fasta):
    """
    Spustí MAFFT pro vícenásobné sekvenční zárovnání.

    Parametry:
    input_fasta (str): Cesta k vstupnímu FASTA souboru.
    output_fasta (str): Cesta k výstupnímu FASTA souboru.
    """
    command = f"mafft --auto {input_fasta} > {output_fasta}"
    subprocess.run(command, shell=True, check=True)

def run_fasttree(input_fasta, output_tree):
    """
    Spustí FastTree pro konstrukci fylogenetického stromu.

    Parametry:
    input_fasta (str): Cesta k vstupnímu FASTA souboru.
    output_tree (str): Cesta k výstupnímu souboru stromu ve formátu Newick.
    """
    command = f"fasttree -nt {input_fasta} > {output_tree}"
    subprocess.run(command, shell=True, check=True)

def plot_tree(tree_file):
    """
    Načte strom z Newick souboru a vykreslí ho.

    Parametry:
    strom_soubor (str): Cesta k Newick souboru.
    """
    tree = Tree(tree_file)

    nstyle = NodeStyle()
    nstyle["shape"] = "circle"
    nstyle["size"] = 5
    nstyle["fgcolor"] = "blue"
    nstyle["vt_line_width"] = 2
    nstyle["hz_line_width"] = 2
    nstyle["vt_line_type"] = 0 
    nstyle["hz_line_type"] = 0

    for node in tree.traverse():
        node.set_style(nstyle)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale = 120 
    ts.branch_vertical_margin = 10  
    ts.title.add_face(TextFace("Core Fylogenetický Strom", fsize=20), column=0)
    tree.show(tree_style=ts)

def main():
    """
    Hlavní funkce pro spuštění celého pipeline na konstrukci fylogenetického stromu core genomů.
    """
    gff_directory = '/cesta/k/souborům/roary/gffs'
    fasta_directory = '/cesta/k/fasta/souborům/'
    core_genes_file = '/cesta/k/souborům/output/core_genes.txt'
    concatenated_fasta = '/cesta/k/souborům/output/concatenated_core_genes.fasta'
    aligned_fasta = '/cesta/k/souborům/output/aligned_core_genes.fasta'
    output_tree = '/cesta/k/souborům/output/core_pangenome_tree.newick'
    tree_file = "/cesta/k/souborům/output/core_pangenome_tree.newick"

    core_genes = load_core_genes(core_genes_file)
    concatenate_core_genes(core_genes, gff_directory, fasta_directory, concatenated_fasta)
    run_mafft(concatenated_fasta, aligned_fasta)
    run_fasttree(aligned_fasta, output_tree)
    plot_tree(tree_file)

    print("Konstrukce Core Fylogenetického Stromu Genomů Dokončena.")

if __name__=="__main__":
    main()
