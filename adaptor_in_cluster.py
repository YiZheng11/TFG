import os
import re

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
            print(line.strip())
            nHits = int(match.group(1))
            break

    return nHits

def save_adaptor_info(file):
    all_adaptor = []
    with open(f"{file}", "r") as adaptor_file:
        adaptor_file.readline()
        for line in adaptor_file:
            record = line.split("\t")[1]
            start = line.split("\t")[2]
            end = line.split("\t")[3]
            strand = line.split("\t")[5]
            hit = (record, start, end, strand)
            all_adaptor.append(hit)
    return all_adaptor

if __name__ == "__main__":
    adaptor_file = os.path.join(os.getcwd(), "adjacent_genes_DUF4123.tsv")
    adaptor_info = save_adaptor_info(adaptor_file)
    all_records = set()
    in_cluster = set()
    for i in range(18):
        n = i + 1
        directory = os.path.join(os.getcwd(), "rizobia_cluster", f"rizobia_cluster_{n}")
        file = os.path.join(directory, "clusterblast_output.txt")
    
        clusters = parse_clusters(file)
        
        # To look at the results
        for adaptor_hit in adaptor_info:
            record = adaptor_hit[0]
            start = adaptor_hit[1]
            end = adaptor_hit[2]
            strand = adaptor_hit[3]
            all_records.add(record)
            for i, cluster in enumerate(clusters):
                index = cluster[2].find(" ")
                cluster_record = (cluster[2]).strip()[index+1:-2]
                if cluster_record in record:
                    print(f"Cluster {n}:")
                    nHits = check_BLAST_hits_in_cluster(cluster)
                    
                    in_cluster.add(record)
                
    not_in_cluster = all_records - in_cluster
    print(len(all_records))
    print(len(in_cluster))
    print(not_in_cluster)
