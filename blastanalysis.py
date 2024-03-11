import os
import subprocess

def blast_data(query, db_name, output_path):
    command = f"tblastn -query {query} -db {db_name} -outfmt 5 -max_target_seqs '2000' -out {output_path}.xml"
    #command = f"tblastn -query {query} -db {db_name} -outfmt '6 qseqid sseqid pident qcovs evalue sseq' -max_target_seqs '2000' -out {output_path}.tsv"
    
    if command is not None:
        subprocess.Popen(command, shell=True).wait()
        print(f"BLAST has been successful. Check out {output_path} where the results are saved.")
    else:
        print("Error, tblastn couldn't be done, please check.")

#def filter_blast(blast_data, output_filtered, minimumIdentity, minimumCoverage, maximumEValue):
#    with open(blast_data, 'r') as input_file, open(output_filtered, 'w') as output_file:
#        for line in input_file:
#            columns = line.split()
#            identity = float(columns[3])
#            coverage = float(columns[4])
#            evalue = float(columns[5])
                
#            if identity >= minimumIdentity and coverage >= minimumCoverage and evalue <= maximumEValue:
#                output_file.write(line)


if __name__ == "__main__":
    query = f"{os.getcwd()}/T6SS_genes_2.fasta"
    db_name = f"{os.getcwd()}/rizobia-data/rizobia_blastdb"
    
    output_dir = "rizobia-blast"
    os.makedirs(output_dir, exist_ok=True)
    output = f"{os.getcwd()}/{output_dir}/tblastnResult_7"
#    output_filtered = "tblastnResult_Filtered_7"
    
    #minimumIdentity = 30
    #minimumCoverage = 50
    #maximumEValue = 1e-5
    
    blast_data(query, db_name, output)
#    filter_blast(output, output_filtered, minimumIdentity, minimumCoverage, maximumEValue)
