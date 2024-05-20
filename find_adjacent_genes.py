import os
import pickle

from Bio import SeqIO, SearchIO
from pathlib import Path

class Accession:
    def __init__(self, accession, path_filename):
        self.accession = accession
        self.path = path_filename
        self.name = ""
        self.loci = []
        self.isplasmid = []
        self.host = ""
        #level of assembly? -> Could not find it in gb file
        
        self.nodA = False
        self.nodB = False
        self.nodC = False
        self.nif = False
        self.tssB = False
        self.tssC = False
        
        self.DUF1795 = False
        self.DUF4123 = False
        self.DUF2169 = False
        self.PRK06147 = False
        self.DUF2875 = False
    
        records = parse(self.path, "genbank")
        for record in records:
            self.loci.append(record.id)
            for feature in record.features:
                if feature.type == "source":
                    plasmid = feature.qualifiers.get("plasmid", None)
                    if plasmid:
                        self.isplasmid.append(True)
                    else:
                        self.isplasmid.append(False)
                    
                    host = feature.qualifiers.get("host", None)
                    if host:
                        self.host = host[0]
                    else:
                        self.host = "Unknown"
                    
                    strain = feature.qualifiers.get("strain", None)
                    if strain:
                        strain = strain[0]
                
                    isolate = feature.qualifiers.get("isolate", None)
                    if isolate:
                        isolate = isolate[0]
        
            species = record.annotations.get("organism")
            if strain and strain not in species:
                species += f" strain {strain}"
            if isolate and isolate not in species:
                species += f" isolate {isolate}"
            self.name = species

def load_pickle_file(pickle_file):
    with open(f"{pickle_file}", 'rb') as file:
        loaded_accessions = pickle.load(file)
    return loaded_accessions

def parse_blast(blast_xml):
    qresults = SearchIO.parse(f"{blast_xml}", "blast-xml")
    ids = []
    for qresult in qresults:
        for hit in qresult.hits:
            for hsp in hit.hsps:
                query_id = qresult.id
                query_start = hsp.query_start + 1
                query_end = hsp.query_end
                ids.append([query_id, query_start, query_end])
    return ids

def parse_genbank(genbank_path):
    gene_positions = []
    for record in SeqIO.parse(genbank_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                GO_function = feature.qualifiers.get("GO_function", [""])[0]
                if GO_function:
                    gene_positions.append((feature.location.start, feature.location.end, GO_function))
    return gene_positions

def find_adjacent_genes(gene_positions, target_start, target_end):
    sorted_genes = sorted(gene_positions, key=lambda x: x[0])
    downstream = []
    upstream = []
    target_index = None

    # Find the index of the target gene in sorted genes
    for i, (start, end, GO) in enumerate(sorted_genes):
        if start <= target_start and end >= target_end:
            target_index = i
            break
    if target_index is not None:
        # Get up to 3 genes to the left
        upstream = sorted_genes[target_index-3:target_index]
        # Get up to 3 genes to the right
        downstream = sorted_genes[target_index+1:target_index+4]

    left_adjacent_genes = [GO for start, end, GO in upstream]
    right_adjacent_genes = [GO for start, end, GO in downstream]

    return left_adjacent_genes, right_adjacent_genes

if __name__ == "__main__":
    rizobia_blast = os.path.join(os.getcwd(), "rizobia-blast")
    DUF1795_file = os.path.join(rizobia_blast, "rpstblastnResult_DUF1795.xml")
    
    pickle_file = os.path.join(os.getcwd(), "accessions_PRUEBA.pkl")
    
    objs = load_pickle_file(pickle_file)
    for obj in objs:
        for query in parse_blast(DUF1795_file):
            for locus_id in obj.loci:
                if query[0] in locus_id:
                    genbank_path = Path(obj.path)
                    target_start = query[1]
                    target_end = query[2]
                    gene_positions = parse_genbank(genbank_path)
                    upstream_genes, downstream_genes = find_adjacent_genes(gene_positions, target_start, target_end)
                    print(genbank_path)
                    print(f"Target gene start: {target_start}")
                    print(f"Target gene end: {target_end}")
                    print(f"Left adjacent genes: {upstream_genes}")
                    print(f"Right adjacent genes: {downstream_genes}")