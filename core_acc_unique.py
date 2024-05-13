import os 

def parse_cdhit_clstr_file(clstr_file):
    clusters = {}
    current_cluster_id = None
    with open(clstr_file, 'r') as file:
        for line in file:
            if line.startswith('>Cluster'):
                current_cluster_id = line.strip().split(' ')[1]
                clusters[current_cluster_id] = []
            else:
                parts = line.strip().split(',')
                seq_id = parts[1].strip().split('...')[0].strip('>')
                clusters[current_cluster_id].append(seq_id)
    return clusters

def categorize_clusters(clusters, total_genomes, core_threshold=0.99):
    core_genome = []
    accessory_genome = []
    unique_genome = []

    for cluster_id, members in clusters.items():
        unique_genomes = set(m.split('_')[0] for m in members) 
        genome_count = len(unique_genomes)
        
        if genome_count == total_genomes:
            core_genome.append(cluster_id)
        elif genome_count == 1:
            unique_genome.append(cluster_id)
        else:
            if genome_count / total_genomes >= core_threshold:
                core_genome.append(cluster_id)
            else:
                accessory_genome.append(cluster_id)

    return core_genome, accessory_genome, unique_genome


def main():
    clstr_dir = '/Users/aijan/Desktop/BP/program/clstr'
    total_genomes = 10  
    all_clusters = {}

    # Aggregate all clusters from .clstr files
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
