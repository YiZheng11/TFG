import os
import pickle
import pdb

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
        #self.seq = ""
        
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

def parse_genbank(genbank_path, locus):
    gene_positions = []
    for record in SeqIO.parse(genbank_path, "genbank"):
        record_id = record.id
        if record_id == locus:
            for feature in record.features:
                if feature.type == "CDS":
                    start = feature.location.start
                    end = feature.location.end
                    strand = feature.location.strand
                    GO_function = feature.qualifiers.get("GO_function", None)
                    if GO_function:
                        gene_positions.append((start, end, GO_function, strand))
                    else:
                        protein = feature.qualifiers.get("product", None)
                        if protein:
                            gene_positions.append((start, end, protein, strand))
                        else:
                            gene_positions.append((start, end, ["Unknown"], strand))
            break
    return gene_positions

def find_adjacent_genes(gene_positions, target_start, target_end):
    sorted_genes = sorted(gene_positions, key=lambda x: x[0])
    downstream = []
    upstream = []
    target = []
    target_index = None
    
    # Find the index of the target gene in sorted genes
    for i, (start, end, GO, strand) in enumerate(sorted_genes):
        if start <= target_start and end >= target_end:
            target_index = i
            break
    if target_index is not None:
        # Get up to 6 genes o the left
        upstream = sorted_genes[target_index-6:target_index]
        # Get up to 6 genes o the right
        downstream = sorted_genes[target_index+1:target_index+7]
        target = sorted_genes[target_index:target_index+1]

    upstream_genes = [(GO, strand) for start, end, GO, strand in upstream]
    downstream_genes = [(GO, strand) for start, end, GO, strand in downstream]
    target_gene = [(GO, strand) for start, end, GO, strand in target]
    
    return upstream_genes, downstream_genes, target_gene

def generate_table(pickle_file, blast_xml, output_name):
    objs = load_pickle_file(pickle_file)
    qresults = parse_blast(blast_xml)
    
    with open(output_name, "w") as table:
        table.write("Organism\tRecord\tStart\tEnd\tFunction\tStrand\tUpstream genes\tDownstream genes\tAccession\n")
        for qresult in qresults:
            for obj in objs:
                for locus in obj.loci:
                    if locus == qresult[0]: #qresult is a list. qresult = [id, start, end]
                        species = obj.name
                        accession = obj.accession
                        genbank_path = Path(obj.path)
                        target_start = qresult[1]
                        target_end = qresult[2]
                        gene_positions = parse_genbank(genbank_path, locus)
                    
                        if gene_positions:
                            upstream, downstream, target = find_adjacent_genes(gene_positions, target_start, target_end)
                            if target:
                                table.write(f"{species}\t{locus}\t{target_start}\t{target_end}\t{target[0][0]}\t{target[0][1]}\t{upstream}\t{downstream}\t{accession}\n")
                            else:
                                table.write(f"{species}\t{locus}\t{target_start}\t{target_end}\tNo info\tNo info\t{upstream}\t{downstream}\t{accession}\n")
                        else:
                            table.write(f"{species}\t{locus}\t{target_start}\t{target_end}\t\t\tNO CDS\tNO CDS\t{accession}\n")
                        break

if __name__ == "__main__":
    rizobia_blast = os.path.join(os.getcwd(), "rizobia-blast")
    DUF1795_file = os.path.join(rizobia_blast, "rpstblastnResult_DUF2169.xml")
    
    pickle_file = os.path.join(os.getcwd(), "accessions.pkl")
    
    generate_table(pickle_file, DUF1795_file, "adjacent_genes_DUF2169.tsv")
