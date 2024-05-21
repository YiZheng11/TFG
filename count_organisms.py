import pdb
import os
import pickle
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

def create_objs(directory, pickle_file):
    "Create objects of class Accession for each GenBank file in <directory> and save them in <pickle_file>"
    directory = Path(directory)
    objs = list()
    for gb_file in Path(directory).rglob("*.gbff"):
        print(f"Creating obj for {gb_file.name}")
        objs.append(Accession(gb_file.name, gb_file))
        
    with open(f"{pickle_file}", 'wb') as file:
        pickle.dump(objs, file)
        
def load_pickle_file(pickle_file):
    "Load <pickle_file>"
    with open(f"{pickle_file}", 'rb') as file:
        loaded_accessions = pickle.load(file)
    return loaded_accessions

def parse_blast(blast_xml, pickle_file):
    "Parse through <blast_xml> and change object attributes saved in <pickle_file>"
    nodA = set()
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
                        obj.nodA = True
                        nodA.add(obj.name)
                        
    with open(f"{pickle_file}", 'wb') as file:
        pickle.dump(objs, file)

    with open("nodA_species", "w") as nodA_file:
        for organism in nodA: #Remember to change
            nodA_file.write(organism+"\n")
            
def parse_rpsblast(blast_xml, pickle_file):
    "Parse through <blast_xml> and change object attributes saved in <pickle_file>"
    DUF2875 = set()
    i=0
    qresults = SearchIO.parse(f"{blast_xml}", "blast-xml")
    for qresult in qresults:
        for obj in objs:
            for locus_id in obj.loci:
                if qresult.id in locus_id:
                    print(qresult.id, locus_id, obj.name)
                    i += 1
                    print(f"nHits: {i}", end="\r")
                    obj.DUF2875 = True
                    DUF2875.add(obj.name)
                        
    with open(f"{pickle_file}", 'wb') as file:
        pickle.dump(objs, file)

    with open("DUF2875_species", "w") as DUF2875_file:
        for organism in DUF2875:
            DUF2875_file.write(organism+"\n")
            
def table_genes(objs):
    "Saves information from <objs> in a txt file and generates a table with the information"
    nodABC = set(); nif = set(); tssBC = set()
    rizobio = set(); rizobio_T6SS = set()
    all_species = set()
    
    DUF1795 = set(); DUF4123 = set(); DUF2169 = set(); PRK06147 = set(); DUF2875 = set()
    
    for obj in objs:
        all_species.add(obj.name)
        
        #If it nodulates
        if obj.nodA == True and obj.nodB == True and obj.nodC == True:
            nodABC.add(obj.name)
        
        #If it's able to fix nitrogen
        if obj.nif == True:
            nif.add(obj.name)
        
        #If it has T6SS
        if obj.tssB == True and obj.tssC == True:
            tssBC.add(obj.name)
        
        #If it nodulates and fix nitrogen = rizobio
        if obj.nodA == True and obj.nodB == True and obj.nodC == True and obj.nif == True:
            rizobio.add(obj.name)
        
        #If it's a rizobio with T6SS
        if obj.nodA == True and obj.nodB == True and obj.nodC == True and obj.nif == True and obj.tssB == True and obj.tssC == True:
            rizobio_T6SS.add(obj.name)
            
            #If the rizobio with T6SS has a certain adaptor
            if obj.DUF1795 == True:
                DUF1795.add(obj.name)
                
            if obj.DUF4123 == True:
                DUF4123.add(obj.name)
            
            if obj.DUF2169 == True:
                DUF2169.add(obj.name)
            
            if obj.PRK06147 == True:
                PRK06147.add(obj.name)
            
            if obj.DUF2875 == True:
                DUF2875.add(obj.name)
            
    #rizobio = nodABC.intersection(nif)
    #rizobio_T6SS = nodABC.intersection(nif, tssBC)
        
    with open("species_nodABC", "w") as nodABC_file:
        nodABC_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in nodABC:
            nodABC_file.write(organism+"\t")
            count = 0 #Keep track of accessions with the same species name
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        nodABC_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        nodABC_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
            
    with open("species_nif", "w") as nif_file:
        nif_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in nif:
            nif_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        nif_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        nif_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
                    
    with open("species_tssBC", "w") as tssBC_file:
        tssBC_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in tssBC:
            tssBC_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        tssBC_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        tssBC_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
                    
    with open("species_rizobio", "w") as rizobio_file:
        rizobio_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in rizobio:
            rizobio_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        rizobio_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        rizobio_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
            
    with open("species_rizobio_T6SS", "w") as rizobio_T6SS_file:
        rizobio_T6SS_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in rizobio_T6SS:
            rizobio_T6SS_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        rizobio_T6SS_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        rizobio_T6SS_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
                    
    with open("species_DUF1795_rizobio_T6SS", "w") as DUF1795_file:
        DUF1795_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in DUF1795:
            DUF1795_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        DUF1795_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        DUF1795_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
                    
    with open("species_DUF4123_rizobio_T6SS", "w") as DUF4123_file:
        DUF4123_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in DUF4123:
            DUF4123_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        DUF4123_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        DUF4123_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
    
    with open("species_DUF2169_rizobio_T6SS", "w") as DUF2169_file:
        DUF2169_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in DUF2169:
            DUF2169_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        DUF2169_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        DUF2169_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
    
    with open("species_PRK06147_rizobio_T6SS", "w") as PRK06147_file:
        PRK06147_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in PRK06147:
            PRK06147_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        PRK06147_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        PRK06147_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
    
    with open("species_DUF2875_rizobio_T6SS", "w") as DUF2875_file:
        DUF2875_file.write("Organism\tnPlasmids\tHost\tAccession\n")
        for organism in DUF2875:
            DUF2875_file.write(organism+"\t")
            count = 0
            for obj in objs:
                if obj.name == organism:
                    if count == 0:
                        DUF2875_file.write(str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    elif count > 0:
                        DUF2875_file.write("\t"+str(sum(obj.isplasmid))+"\t"+obj.host+"\t"+obj.accession+"\n")
                    count += 1
                    
    with open("species_all", "w") as all_file:
        for organism in all_species:
            all_file.write(organism+"\n")
    
    
    files = ["species_all", "species_nodABC", "species_nif", "species_rizobio", "species_tssBC", "species_rizobio_T6SS"]
    files_adaptor = ["species_rizobio_T6SS", "species_DUF1795_rizobio_T6SS", "species_DUF4123_rizobio_T6SS",
                     "species_DUF2169_rizobio_T6SS", "species_PRK06147_rizobio_T6SS", "species_DUF2875_rizobio_T6SS"]
    table1 = []
    table2 = []
    
    headers1 = ['Genus', 'Total', 'nodABC', 'nifH', 'rizobios', 'tssBC', 'rizobios con T6SS']
    headers2 = ['Genus', 'rizobios con T6SS', 'DUF1795', 'DUF4123', 'DUF2169', 'PRK06147', 'DUF2875']
    index = ["Rhizobium", "Bradyrhizobium", "Mesorhizobium", "Sinorhizobium"]
    for file in files:
        column = []
        species_file = os.path.join(os.getcwd(), f"{file}")
        rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium, total, extra = number_genus(species_file)
        column = [rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium]
        table1.append(column)
    df1 = pd.DataFrame(zip(*table1))
    df1.insert(0, "Genus", index)
    df1.to_csv("table1.tsv", sep = '\t', header=headers1, index=False)
    
    for file in files_adaptor:
        column = []
        species_file = os.path.join(os.getcwd(), f"{file}")
        rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium, total, extra = number_genus(species_file)
        column = [rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium]
        table2.append(column)
    df2 = pd.DataFrame(zip(*table2))
    df2.insert(0, "Genus", index)
    df2.to_csv("table2.tsv", sep = '\t', header=headers2, index=False)
        
def number_genus(file):
    "Counts species name stored in <file> for each genera"
    rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium = 0, 0, 0, 0
    total, extra = 0, 0
    with open(f"{file}", "r") as species_file:
        for line in species_file:
            total += 1
            column_1 = line.split("\t")[0]
            if "Rhizobium" in column_1 or "Agrobacterium" in column_1:
                rhizobium += 1
            elif "Bradyrhizobium" in column_1:
                bradyrhizobium += 1
            elif "Mesorhizobium" in column_1:
                mesorhizobium += 1
            elif "Sinorhizobium" in column_1 or "Ensifer" in column_1:
                sinorhizobium += 1
            else:
                extra += 1
    return rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium, total, extra

if __name__ == "__main__":
    rizobia_blast = os.path.join(os.getcwd(), "rizobia-blast")
    nodA_file = os.path.join(rizobia_blast, "tblastnResult_nodA.xml")
    nodB_file = os.path.join(rizobia_blast, "tblastnResult_nodB.xml")
    nodC_file = os.path.join(rizobia_blast, "tblastnResult_nodC.xml")
    nifH_file = os.path.join(rizobia_blast, "tblastnResult_nifH.xml")
    tssB_file = os.path.join(rizobia_blast, "tblastnResult_tssB.xml")
    tssC_file = os.path.join(rizobia_blast, "tblastnResult_tssC.xml")
    
    DUF1795_file = os.path.join(rizobia_blast, "rpstblastnResult_DUF1795.xml")
    DUF4123_file = os.path.join(rizobia_blast, "rpstblastnResult_DUF4123.xml")
    DUF2169_file = os.path.join(rizobia_blast, "rpstblastnResult_DUF2169.xml")
    PRK06147_file = os.path.join(rizobia_blast, "rpstblastnResult_PRK06147.xml")
    DUF2875_file = os.path.join(rizobia_blast, "rpstblastnResult_DUF2875.xml")
    
    data_dir = os.path.join(os.getcwd(), "rizobia-data")
    pickle_file = os.path.join(os.getcwd(), "accessions_PRUEBA.pkl")
    
    #create_objs(data_dir, pickle_file)
    objs = load_pickle_file(pickle_file)
    #parse_rpsblast(DUF2875_file, pickle_file) #Remember to change
    table_genes(objs)
