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
        unique_genomes = set(m.split('_')[0] for m in members)  # Assuming genome ID can be derived from sequence IDs
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

def main(clstr_file, total_genomes):
    clusters = parse_cdhit_clstr_file(clstr_file)
    core, accessory, unique = categorize_clusters(clusters, total_genomes)
    
    print(f"Core Genome Clusters: {len(core)}")
    print(f"Accessory Genome Clusters: {len(accessory)}")
    print(f"Unique Genome Clusters: {len(unique)}")

clstr_file = '/Users/aijan/Desktop/MAIN PIPELINE/output/cdhit_output.clstr'
total_genomes = 10
main(clstr_file, total_genomes)