import os
import subprocess

def blast_data(query, db_name, output_file):
    command = f"tblastn -query {query} -db {db_name} -outfmt '6 qseqid sseqid sseq pident qcovs evalue ' -max_target_seqs '2000' -out {output_file}"
    
    if command is not None:
        subprocess.Popen(command, shell=True).wait()
        print(f"BLAST has been successful. Check the file {output_file} where the results are shown.")
    else:
        print("Error, tblastn couldn't be done, please check.")

def filter_blast(blast_data, output_filtered, minimumIdentity, minimumCoverage, maximumEValue):
    with open(blast_data, 'r') as input_file, open(output_filtered, 'w') as output_file:
        for line in input_file:
            columns = line.split()
            identity = float(columns[3])
            coverage = float(columns[4])
            evalue = float(columns[5])
                
            if identity >= minimumIdentity and coverage >= minimumCoverage and evalue <= maximumEValue:
                output_file.write(line)
        

if __name__ == "__main__":
    query = f"{os.getcwd()}/T6SS_genes_2.fasta"
    db_name = f"{os.getcwd()}/data/bradyrhizobia_blastdb"
    output = "tblastnResult_6.tsv"
    output_filtered = "tblastnResult_Filtered_6.tsv"
    
    minimumIdentity = 30
    minimumCoverage = 50
    maximumEValue = 1e-5
    
    blast_data(query, db_name, output)
    filter_blast(output, output_filtered, minimumIdentity, minimumCoverage, maximumEValue)