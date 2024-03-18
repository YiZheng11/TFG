import os
import subprocess

def blast_data(query, db_name, output_path):
    command = f"tblastn -query {query} -db {db_name} -outfmt 5 -max_target_seqs '2000' -out {output_path}"
    
    if command is not None:
        subprocess.Popen(command, shell=True).wait()
        print(f"BLAST has been successful. Check out {output_path} where the results are saved.")
    else:
        print("Error, tblastn couldn't be done, please check.")

def filter_blast(blast_xml, output_file, minimumIdentity, minimumCoverage, maximumEValue):
    with open(output_file, 'w') as out_handle:
        out_handle.write("query_id\thit_id\tpident\tqcovs\tevalue\n")
        result_handle = open(blast_xml)
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= maximumEValue and hsp.align_length / blast_record.query_length * 100 >= minimumCoverage and hsp.identities / hsp.align_length * 100 >= minimumIdentity:
                        out_handle.write(f"{blast_record.query_def}\t{alignment.hit_id}\t{hsp.identities / hsp.align_length * 100:.2f}\t{hsp.align_length / blast_record.query_length * 100:.2f}\t{hsp.expect}\n")
        result_handle.close()

if __name__ == "__main__":
    query = f"{os.getcwd()}/T6SS_genes_2.fasta"
    db_name = f"{os.getcwd()}/rizobia-data/rizobia_blastdb"
    
    output_dir = "rizobia-blast"
    os.makedirs(output_dir, exist_ok=True)
    output = f"{os.getcwd()}/{output_dir}/tblastnResult_T6SS_I.xml"
    output_filtered = "tblastnResult_T6SS_I_Filtered"
    
    #minimumIdentity = 30
    #minimumCoverage = 50
    #maximumEValue = 1e-5
    
    blast_data(query, db_name, output)
    filter_blast(output, output_filtered, minimumIdentity, minimumCoverage, maximumEValue)
