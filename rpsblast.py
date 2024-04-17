import os
import subprocess
import tempfile

from Bio import SeqIO
from pathlib import Path

def profile_db(PSSM_file, output):
    "Create profile db with the PSSMs named in <PSSM_file> and saved it in <output>"
    command = f"makeprofiledb -in {PSSM_file} -out {output}"

    if command is not None:
        subprocess.Popen(command, shell=True).wait()
        print(f"Check out {output} where the db is saved.")
    else:
        print("Error, makeprofiledb couldn't be done, please check.")

def create_multifasta(directory):
    "Create a temporary multifasta file from all translated sequences in the directory."
    #temp_multifasta = tempfile.NamedTemporaryFile(delete=False)
    directory = Path(directory) # QUITAR LUEGO
    gb_files = directory.rglob("*.gbff")

    with open("temp_multifasta2", "w") as multifasta_file:
        for gb_file in gb_files:
            for record in SeqIO.parse(gb_file, "genbank"):
                locus_id = record.id
                for feature in record.features:
                    if feature.type == "CDS":
                        protein_id = feature.qualifiers.get("protein_id", ["Unknown"])[0]
                        translation = feature.qualifiers.get("translation", [""])[0]
                        if translation:
                            multifasta_file.write(f">{locus_id} {protein_id}\n{translation}\n")
                        else:
                            try:
                                location = feature.location
                                nucleotide_seq = location.extract(record).seq
                                translation = nucleotide_seq.translate()
                                multifasta_file.write(f">{locus_id} {protein_id}\n{translation}\n")
                            except:
                                multifasta_file.write(f">{locus_id} {protein_id}\n{translation}\n")

    return None

def rpsblast(directory, profile_db, output_name):
    "Perform RPS BLAST using sequences extracted from GenBank files"
    directory = Path(directory)
    #temp_multifasta_file = create_multifasta(directory)
    temp_multifasta_file = f"{os.getcwd()}/temp_multifasta2"

    command = ["rpsblast", "-query", temp_multifasta_file, "-db", profile_db, "-out", output_name, "-outfmt", "5", "-max_target_seqs", "8000", "-evalue", "0.01"]
    #command = f"rpsblast -query {temp_multifasta_file} -db {profile_db} -out {output_name} -outfmt 5 -max_target_seqs '8000' -evalue 0.01"
    
    with open(temp_multifasta_file, 'rb') as input_file:
        process = subprocess.Popen(command, stdin=input_file, stderr=subprocess.PIPE)
        _, stderr = process.communicate()

    # Remove the temporary multifasta file
    #Path("tmpebuzvxfx").unlink()
    return None

if __name__ == "__main__":
    PSSM_file = os.path.join(os.getcwd(), "DUF1795.pn")
    output_dir_profile = os.path.join(os.getcwd(), "rizobia-profiledb")
    os.makedirs(output_dir_profile, exist_ok=True)
    output_name_profile = "DUF1795_profiledb"

    dir_protein = os.path.join(os.getcwd(), "rizobia-data")
    output_rps = os.path.join(os.getcwd(), "rizobia-blast", "rpsblastResult_DUF1795.xml")

    profile_db(PSSM_file, os.path.join(output_dir_profile, output_name_profile))
    rpsblast(dir_protein, f"{os.getcwd()}/rizobia-profiledb/DUF1795_profiledb", output_rps)