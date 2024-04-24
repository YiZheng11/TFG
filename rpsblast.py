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

def cds_feature_to_record(cds):
    seq = Seq(cds.qualifiers.get("translation", [""])[0])
    protein_id = cds.qualifiers.get("protein_id", [""])[0]
    locus_tag = cds.qualifiers.get("locus_tag", [""])[0]
    return SeqRecord(seq, id=locus_tag, description=protein_id)

def fasta_generator(directory, genus):
    for gb_file in Path(directory).rglob("*.gbff"):
        print(f"For {gb_file.name}")
        for record in parse(gb_file, "genbank"):
            multifasta = ""
            for cds in filter(lambda x: x.type == "CDS", record.features):
                protein = cds_feature_to_record(cds)
                if len(protein.seq) > 0:
                    multifasta += protein.format("fasta")
            yield multifasta
    return None

def rpsblast(directory, profile_db, output_name, genus):
    "Perform RPS BLAST using sequences extracted from GenBank files"
    directory = Path(directory)
    # temp_multifasta_file = create_multifasta(directory, genus)

    command = ["rpsblast+"]
    command += ["-query",  "-"]
    command += ["-db",  profile_db]
    command += ["-outfmt", "5"]
    command += ["-evalue",  "0.001"]
    command += ["-max_hsps",  "1"]

    queries = 0
    hits = 0
    qresults = []
    for fasta in fasta_generator(directory, genus):
        if fasta:
            proc = subprocess.Popen(command, text=True, stdin=PIPE, stdout=PIPE)
            xmlout, _ = proc.communicate(input=fasta)
            for qresult in SearchIO.parse(io.StringIO(xmlout), "blast-xml"):
                queries += 1
                if len(qresult.hits) > 0:
                    hits += len(qresult.hits)
                    qresults += [qresult]
        print(f"nQueries: {queries} nHits: {hits}", end="\r")
    SearchIO.write(qresults, "myresults.xml", "blast-xml")
    return None

if __name__ == "__main__":
    PSSM_file = os.path.join(os.getcwd(), "DUF1795.pn")
    output_dir_profile = os.path.join(os.getcwd(), "rhizobia-profiledb")
    os.makedirs(output_dir_profile, exist_ok=True)
    output_name_profile = "DUF1795_profiledb"

    dir_protein = "/home/lewis/ISLU/data/Bradyrhizobium-data/"
    output_rps = os.path.join(os.getcwd(), "rpsblastResult_DUF1795_Mesorhizobium.xml")

    # profile_db(PSSM_file, os.path.join(output_dir_profile, output_name_profile))
    rpsblast(dir_protein, "/home/lewis/Data/cdd/duf1795", output_rps, "Bradyrhizobium")
