import os
import io
import subprocess
from subprocess import PIPE

from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from Bio import SearchIO

def profile_db(PSSM_file, output):
    "Create profile db with the PSSMs named in <PSSM_file> and saved it in <output>"
    command = f"makeprofiledb -in {PSSM_file} -out {output}"
    subprocess.Popen(command, shell=True).wait()
    print(f"Check out {output} where the db is saved.")
    
def fasta_generator(directory):
    "Iterates over FNA files in <directory> and yield it as a string"
    #for fasta_file in Path(directory).rglob("*.fna"):
    #print(f"For {fasta_file.name}")
    with open(directory, "r") as f:
        fasta_content = f.read()
    yield fasta_content
    return None

def rpstblastn(directory, profile_db, output_rps):
    "Perform RPS TBLASTN using sequences extracted from FNA files"
    #directory = Path(directory)
    
    command = ["rpstblastn"]
    command += ["-query", "-"]
    command += ["-db", profile_db]
    command += ["-outfmt", "5"]
    #command += ["max_target_seqs", "8000"]
    command += ["-evalue", "0.00001"]
    
    queries = 0
    hits = 0
    qresults = []
    for fasta in fasta_generator(directory):
        if fasta:
            process = subprocess.Popen(command, text=True, stdin=PIPE, stdout=PIPE)
            xmlout, _ = process.communicate(input=fasta)
            for qresult in SearchIO.parse(io.StringIO(xmlout), "blast-xml"):
                queries += 1
                
                if len(qresult.hits) > 0:
                    hits += len(qresult.hits)
                    qresults += [qresult]
        #print(f"nQueries: {queries} nHits: {hits}", end="\r")
        print("nQueries:", queries, "nHits:", hits)
    SearchIO.write(qresults, output_rps, "blast-xml")
    return None

if __name__ == "__main__":
    PSSM_file = os.path.join(os.getcwd(), "DUF1795.pn")
    output_dir_profile = os.path.join(os.getcwd(), "rizobia-profiledb")
    os.makedirs(output_dir_profile, exist_ok=True)
    output_name_profile = "DUF1795_profiledb"

    data_dir = os.path.join(os.getcwd(), "rizobia-data")
    output_rps = os.path.join(os.getcwd(), "rizobia-blast", "rpstblastnResult_DUF1795_PRUEBA.xml")
    
    mesorhizobium = os.path.join(os.getcwd(), "rizobia-data", "GCF_022347545.1_ASM2234754v1_genomic.fna")

    #profile_db(PSSM_file, os.path.join(output_dir_profile, output_name_profile))
    rpstblastn(mesorhizobium, os.path.join(output_dir_profile, output_name_profile), output_rps)
    
