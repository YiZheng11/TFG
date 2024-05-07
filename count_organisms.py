import os
import pickle

from pathlib import Path
from Bio.Blast import NCBIXML
from Bio.SeqIO import parse
from Bio import SeqIO

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

if __name__ == "__main__":
    nodABC_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_nodABC.xml"
    nifH_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_nifH.xml"
    T6SS_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_T6SS.xml"

    data_dir = os.path.join(os.getcwd(), "rizobia-data")
    
    genera = ["Rhizobium", "Bradyrhizobium", "Mesorhizobium", "Sinorhizobium"]
    create_objs(data_dir)
