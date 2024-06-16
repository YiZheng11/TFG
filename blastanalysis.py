import os
import subprocess
import pdb
import time

from Bio.Blast import NCBIXML
from Bio import SeqIO
from urllib.error import HTTPError 

def blast_data(query, db_name, output_path):
    command = f"tblastn -query {query} -db {db_name} -outfmt 5 -max_target_seqs '8000' -out {output_path} -evalue 0.00001"

    try:
        subprocess.Popen(command, shell=True).wait()
        print(f"BLAST has been successful. Check out {output_path} where the results are saved.")
    except:
        print("Error, tblastn couldn't be done, please check.")

def filter_blast(blast_xml, output_file):
    with open(output_file, 'w') as out_handle:
        out_handle.write("query_id\thit_id\thit_def\n")
        result_handle = open(blast_xml)
        blast_record = NCBIXML.parse(result_handle)
        for query in blast_record:
            for alignment in query.alignments:
                for hsp in alignment.hsps:
                    #if (hsp.expect <= maximumEValue) and ((hsp.query_end - hsp.query_start + 1)/ query.query_length * 100 >= minimumCoverage) and (hsp.identities / hsp.align_length * 100 >= minimumIdentity):
                        #annotation_status = check_annotation(f"{os.getcwd()}/bradyrhizobia-data", f"{alignment.hit_id}", f"{hsp.sbjct_start}", f"{hsp.sbjct_end}")
                        #out_handle.write(f"{query.query}\t{alignment.hit_id}\t{hsp.identities / hsp.align_length * 100:.2f}\t{hsp.align_length / query.query_length * 100:.2f}\t{hsp.expect}\t{annotation_status}\n")
                        #out_handle.write(f"{query.query}\t{alignment.hit_id}\t{alignment.hit_def}\t{hsp.identities / hsp.align_length * 100:.2f}\t{(hsp.query_end - hsp.query_start + 1)/ query.query_length * 100:.2f}\t{hsp.expect}\n")
                    out_handle.write(f"{query.query}\t{alignment.hit_id}\t{alignment.hit_def}\n")
        result_handle.close()
    print(f"Filtering done. Check out {output_file} where the filtered results are saved")

def check_annotation(db_directory, hit_id, hit_start, hit_end):
    hit_id = hit_id.split("|")[1] #Example: NZ_CP029603.1
    locus_id = hit_id.split(".")[0] #Example: NZ_CP029603
    
    annotation_status = ""
    for gb_file in os.listdir(db_directory):
        if gb_file.endswith(".gbff"):
            for record in SeqIO.parse(gb_file, "genbank"):
                if record.id == locus_id or record.id == hit_id:
                    for feature in record.features:
                        if feature.type == "CDS":
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            if start <= int(hit_start) <= end or start <= int(hit_end) <= end:
                                # Check if the CDS partially overlaps with the hit range
                                if (start <= int(hit_start) <= end) and (start <= int(hit_end) <= end):
                                    # If partially overlaps, return "Ambiguous"
                                    annotation_status = "Annotated"
                                    print(annotation_status)
                                    return annotation_status
                                else:
                                    # If fully overlaps, return "Annotated"
                                    annotation_status = "Ambiguous"
                                    print(annotation_status)
                                    return annotation_status
    # If no overlapping CDS is found, return "Not Annotated"
    return "Not Annotated"

def filter_blast_with_retry(blast_xml, output_file, minimumIdentity, minimumCoverage, maximumEValue, max_retries=20):
    retries = 0
    success = False

    while not success and retries < max_retries:
        try:
            filter_blast(blast_xml, output_file, minimumIdentity, minimumCoverage, maximumEValue)
            success = True
        except HTTPError as e:
            if e.code == 502:  # HTTPError: Bad Gateway
                print(f"Internal Server Error ({e.code}), retrying again in 60 seconds... Retry attempt: {retries + 1}/{max_retries}")
                time.sleep(60)
                retries += 1
            else:
                raise  # If it is not HTTPError: Bad Gateway, it shows the exception/error captured

    if not success:
        print(f"Failed to filter the results of {blast_xml} after {max_retries} retries.")


if __name__ == "__main__":
    gene = "tssC"
    
    query = f"{os.getcwd()}/{gene}.fasta"
    db_name = f"{os.getcwd()}/rizobia-data/rizobia_blastdb"
    
    output_dir = "rizobia-blast"
    os.makedirs(output_dir, exist_ok=True)
    output = f"{os.getcwd()}/{output_dir}/tblastnResult_{gene}.xml"
    
    output_filtered = f"{os.getcwd()}/{output_dir}/tblastnResult_{gene}_Filtered_noAnnotation"
    
    #minimumIdentity = 30
    #minimumCoverage = 50
    #maximumEValue = 1e-5
    
    blast_data(query, db_name, output)
    #filter_blast(output, output_filtered)    
