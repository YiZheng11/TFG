import os
from Bio.Blast import NCBIXML
from Bio import SeqIO, Entrez

Entrez.email = "yi.zheng@alumnos.upm.es" # your email please
Entrez.api_key = "1c105008de567a0fdcc74eedb9584b2ec109"  # your API key
    
class Organism:
    def __init__(self):
        self.organisms = {}

    def add_organism(self, name, queries=None):
        if queries is None:
            queries = set()  
        if name in self.organisms:
            self.organisms[name].update(queries)
        else:
            self.organisms[name] = set(queries)

    def count_organisms(self, genus, required_queries=None):
        count = 0
        for organism, queries in self.organisms.items():
            if genus in organism and (required_queries is None or queries.intersection(required_queries) == required_queries):
                count += 1
        return count

def parse_xml(filename, organism_collection, queries=None):
    result_handle = open(filename)
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hit_accession = alignment.accession
                species, strain, isolate = get_hit_info(hit_accession)
                organism_name = species
                if strain and strain not in species:
                    organism_name += f" {strain}"
                if isolate and isolate not in species:
                    organism_name += f" {isolate}"
                print(organism_name)
                
                query_name = None
                if queries:
                    for query in queries:
                        if query in blast_record.query:
                            query_name = query
                            organism_collection.add_organism(organism_name, query_name)
                else:
                    organism_collection.add_organism(organism_name)
                
                break  # Solo necesitamos la información del primer HSP
    result_handle.close()

def get_hit_info(hit_accession):
    species = None; strain = None; isolate = None
    try:
        handle = Entrez.efetch(db="nuccore", id=hit_accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        species = record.annotations.get("organism", "Unknown species")
        for feature in record.features:
            if feature.type == "source":
                strain = feature.qualifiers.get("strain", None)
                if strain:
                    strain = strain[0]
                isolate = feature.qualifiers.get("isolate", None)
                if isolate:
                    isolate = isolate[0]
                break
    except Exception as e:
        print("Error:", e)
    return species, strain, isolate

if __name__ == "__main__":
    nodABC_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_nodABC.xml"
    nifH_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_nifH.xml"
    T6SS_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_T6SS.xml"
    
    organism_collection_nifH = Organism()
    parse_xml(nifH_file, organism_collection_nifH)
    genera = ["Rhizobium", "Bradyrhizobium", "Mesorhizobium", "Sinorhizobium"]
    for genus in genera: 
        #organism_collection_nodABC = Organism()
        #parse_xml(nodABC_file, organism_collection_nodABC, data_dir, ["NodA", "NodB", "NodC"])
        #num_organisms_nodABC = organism_collection_nodABC.count_organisms(genus, {"NodA", "NodB", "NodC"})

        num_organisms_nifH = organism_collection_nifH.count_organisms(genus)

        #organism_collection_T6SS = Organism()
        #parse_xml(T6SS_file, organism_collection_T6SS, data_dir)
        #num_organisms_T6SS = organism_collection_T6SS.count_organisms(genus)

        #print(f"{genus} nodABC: {num_organisms_nodABC}")
        print(f"{genus} nifH: {num_organisms_nifH}")
        #print(f"{genus} T6SS: {num_organisms_T6SS}")
