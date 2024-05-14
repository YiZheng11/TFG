import os
import pickle
import pdb
import csv
import pandas as pd

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
        #level of assembly? -> Could not find it in gb file

        self.nodA = False
        self.nodB = False
        self.nodC = False
        self.nif = False
        self.tssB = False
        self.tssC = False
    
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

def table_genes(objs):
    nodABC = set()
    nif = set()
    tssBC = set()
    rizobio = set()
    rizobio_T6SS = set()
    all_species = set()
    
    for obj in objs:
        all_species.add(obj.name)
        if obj.nodA == True and obj.nodB == True and obj.nodC == True:
            nodABC.add(obj.name)
            
        if obj.nif == True:
            nif.add(obj.name)
        
        if obj.tssB == True and obj.tssC == True:
            tssBC.add(obj.name)
            
        if obj.nodA == True and obj.nodB == True and obj.nodC == True and obj.nif == True:
            rizobio.add(obj.name)
        
        if obj.nodA == True and obj.nodB == True and obj.nodC == True and obj.nif == True and obj.tssB == True and obj.tssC == True:
            rizobio_T6SS.add(obj.name)
            
    #rizobio = nodABC.intersection(nif)
    #rizobio_T6SS = nodABC.intersection(nif, tssBC)
        
    with open("species_nodABC", "w") as nodABC_file:
        nodABC_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in nodABC:
            nodABC_file.write(organism+"\t")
            for obj in objs:
                if obj.name == organism:
                    nodABC_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
            
    with open("species_nif", "w") as nif_file:
        nif_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in nif:
            nif_file.write(organism+"\t")
            for obj in objs:
                if obj.name == organism:
                    nif_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    
    with open("species_tssBC", "w") as tssBC_file:
        tssBC_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in tssBC:
            tssBC_file.write(organism+"\t")
            for obj in objs:
                if obj.name == organism:
                    tssBC_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    
    with open("species_rizobio", "w") as rizobio_file:
        rizobio_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in rizobio:
            rizobio_file.write(organism+"\t")
            for obj in objs:
                if obj.name == organism:
                    rizobio_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
            
    with open("species_rizobio_T6SS", "w") as rizobio_T6SS_file:
        rizobio_T6SS_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in rizobio_T6SS:
            rizobio_T6SS_file.write(organism+"\t")
            for obj in objs:
                if obj.name == organism:
                    rizobio_T6SS_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    
    with open("species_all", "w") as all_file:
        for organism in all_species:
            all_file.write(organism+"\n")
    
    
    files = ["species_nodABC", "species_nif", "species_rizobio", "species_tssBC", "species_rizobio_T6SS"]
    table = []
    
    headers = ['Genus', 'nodABC', 'nifH', 'rizobios', 'tssBC', 'rizobios con T6SS']
    index = ["Rhizobium", "Bradyrhizobium", "Mesorhizobium", "Sinorhizobium"]
    for file in files:
        column = []
        species_file = os.path.join(os.getcwd(), f"{file}")
        rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium, total, extra = number_genus(species_file)
        print(total, extra)
        column = [rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium]
        table.append(column)
    df = pd.DataFrame(zip(*table))
    df.insert(0, "Genus", index)
    df.to_csv("table.tsv", sep = '\t', header=headers, index=False)

def number_genus(file):
    rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium = 0, 0, 0, 0
    total, extra = 0, 0
    with open(f"{file}", "r") as species_file:
        for line in species_file:
            total += 1
            column_1 = line.split("\t")[0] #The first column contains the species name
            if "Rhizobium" in column_1 or "Agrobacterium" in column_1:
                rhizobium += 1
            elif "Bradyrhizobium" in column_1:
                bradyrhizobium += 1
            elif "Mesorhizobium" in column_1:
                mesorhizobium += 1
            elif "Sinorhizobium" in column_1 or "Ensifer" in column_1:
                sinorhizobium += 1
            else:
                extra += 1 #Check the lines in the file where the fist column doesn't have the species name
    return rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium, total, extra

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
