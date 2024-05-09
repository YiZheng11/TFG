import os
import pickle

from pathlib import Path
from Bio.SeqIO import parse
from Bio import SeqIO, SearchIO

class Accession:
    def __init__(self, accession, path_filename):
        self.accession = accession
        self.path = path_filename
        self.name = ""
        self.loci = []
        self.isplasmid = []
        self.host = ""
    
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

def create_objs(directory):
    directory = Path(directory)
    objs = list()
    for gb_file in Path(directory).rglob("*.gbff"):
        print(f"Creating obj for {gb_file.name}")
        objs.append(Accession(gb_file.name, gb_file))
        
    with open('accessions.pkl', 'wb') as file:
            pickle.dump(objs, file)

def load_pickle_file(pickle_file):
    with open(f"{pickle_file}", 'rb') as file:
        loaded_accessions = pickle.load(file)
    return loaded_accessions

def parse_blast(blast_xml, pickle_file):
    nodC = set()
    i=0
    qresults = SearchIO.parse(f"{blast_xml}", "blast-xml")
    for qresult in qresults:
        for hit in qresult.hits:
            for obj in objs:
                for locus_id in obj.loci:
                    if hit.accession in locus_id:
                        print(hit.accession, locus_id, obj.name)
                        i += 1
                        print(f"nHits: {i}", end="\r")
                        obj.nodC = True
                        nodC.add(obj.name)
                        
    with open(f"{pickle_file}", 'wb') as file:
        pickle.dump(objs, file)
    
    with open("nodC_species", "w") as nodC_file:
        for organism in nodC:
            nodC_file.write(organism+"\n")

if __name__ == "__main__":
    rizobia_blast = os.path.join(os.getcwd(), "rizobia-blast")
    nodA_file = os.path.join(rizobia_blast, "tblastnResult_nodA.xml")
    nodB_file = os.path.join(rizobia_blast, "tblastnResult_nodB.xml")
    nodC_file = os.path.join(rizobia_blast, "tblastnResult_nodC.xml")
    nifH_file = os.path.join(rizobia_blast, "tblastnResult_nifH.xml")
    T6SS_file = os.path.join(rizobia_blast, "tblastnResult_T6SS.xml")

    data_dir = os.path.join(os.getcwd(), "rizobia-data")
    pickle_file = os.path.join(os.getcwd(), "accessions.pkl")
    
    genera = ["Rhizobium", "Bradyrhizobium", "Mesorhizobium", "Sinorhizobium"]
    #create_objs(data_dir, pickle_file)
    objs = load_pickle_file(pickle_file)
    parse_blast(nodC_file, pickle_file)
