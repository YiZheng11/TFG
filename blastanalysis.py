import os
import subprocess

def blast_data(query, db_name, output_file):
    command = f"tblastn -query {query} -db {db_name} -outfmt '6 qseqid sseq pident qcovs evalue sseqid' -max_target_seqs '2000' -out {output_file}"
    
    if command is not None:
        subprocess.Popen(command, shell=True).wait()
        print(f"BLAST has been successful. Check the {output_file} where the results are shown")
    else:
        print("Error, tblastn couldn't be done, please check.")

if __name__ == "__main__":
    query = f"{os.getcwd()}/T6SS_genes.fasta"
    db_name = f"{os.getcwd()}/data/bradyrhizobia_blastdb"
    output = "tblastnResult.tsv"
    blast_data(query, db_name, output)