import os
import re
import pdb

def parse_clusters(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    clusters = []
    checking_cluster = []

    for line in lines:
        if line.startswith('>>'):
            if checking_cluster:
                clusters.append(checking_cluster)
                checking_cluster = []
            checking_cluster.append(line)
        else:
            if checking_cluster:
                checking_cluster.append(line)

    clusters.append(checking_cluster)
    return clusters

def check_BLAST_hits_in_cluster(cluster):
    nHits = 0
    pattern = re.compile(r"Number of proteins with BLAST hits to this cluster:\s*(\d+)")
    
    for line in cluster:
        match = pattern.search(line)
        if match:
            nHits = int(match.group(1))
            break
    return nHits

def check_location(cluster, adaptor_start, adaptor_end):
    cluster_start = None
    cluster_end = None
    cluster_lines = []
    found = False
    
    for line in cluster:
        if not found:
            if line.startswith("Table of genes, locations, strands and annotations of subject cluster:"):
                found = True
                continue
        
        if found:
            if line.startswith("\n"):
                break
            cluster_lines.append(line)
            
    cluster_start = int(cluster_lines[0].split("\t")[1])
    cluster_end = int(cluster_lines[-1].split("\t")[2])
    
    if cluster_start <= adaptor_start <= cluster_end and cluster_start <= adaptor_end <= cluster_end:
        return True
    else:
        return False

def save_adaptor_info(file):
    all_adaptor = []
    with open(file, "r") as adaptor_file:
        adaptor_file.readline()  # Skip the header line
        for line in adaptor_file:
            adaptor = line.strip().split("\t")
            all_adaptor.append(adaptor)
    return all_adaptor

def main(adaptor_file, output_file):
    adaptor_info = save_adaptor_info(adaptor_file)
    all_records = []
    not_in_cluster = []
    in_cluster = []
    
    for adaptor_hit in adaptor_info:
        record = adaptor_hit[1]
        start = int(adaptor_hit[2])
        end = int(adaptor_hit[3])
        strand = adaptor_hit[5]
        all_records.append(record)
        
        found_in_cluster = False
        n = 0
        for i in range(18):
            n = i + 1
            directory = os.path.join(os.getcwd(), "rizobia_cluster", f"rizobia_cluster_{n}")
            file = os.path.join(directory, "clusterblast_output.txt")
    
            clusters = parse_clusters(file)
            for cluster in clusters:
                index = cluster[2].find(" ")
                cluster_record = cluster[2].strip()[index+1:-2]
                if cluster_record in record:
                    nHits = check_BLAST_hits_in_cluster(cluster)
                    if nHits > 4:
                        adaptor_in = check_location(cluster, start, end)
                        if adaptor_in:
                            found_in_cluster = True
                            in_cluster.append((record, start, end, strand))
                            break
            if found_in_cluster:
                break
        
        if not found_in_cluster:
            not_in_cluster.append(adaptor_hit)
                        
            
    with open(output_file, 'w') as outfile:
        for adaptor in not_in_cluster:
            outfile.write("\t".join(map(str, adaptor)) + "\n")
    
    print(f"Total number of adaptors found initially: {len(all_records)}")
    print(f"Number of those that are in a cluster: {len(in_cluster)}")
    print(f"Number of those that are NOT in a cluster: {len(not_in_cluster)}")
    print(f"Results for adaptors not found in any cluster saved in: {output_file}")

if __name__ == "__main__":
    DUF = "DUF2169"
    output_dir = os.path.join(os.getcwd(), "adaptor-not-in-cluster", DUF)
    os.makedirs(f"{output_dir}", exist_ok=True)
    
    output_file = os.path.join(output_dir, f"adaptor_not_in_cluster.tsv")
    adaptor_file = os.path.join(os.getcwd(), "adjacent-genes", DUF, f"adjacent_genes.tsv")
    
    main(adaptor_file, output_file)
    
