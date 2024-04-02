import os
import subprocess
import pdb
import time

from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
from urllib.error import HTTPError 

Entrez.email = "yi.zheng@alumnos.upm.es" # your email please
Entrez.api_key = "1c105008de567a0fdcc74eedb9584b2ec109"  # your API key

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
        blast_record = NCBIXML.parse(result_handle)
        for query in blast_record:
            for alignment in query.alignments:
                for hsp in alignment.hsps:
                    if (hsp.expect <= maximumEValue) and ((hsp.query_end - hsp.query_start + 1)/ query.query_length * 100 >= minimumCoverage) and (hsp.identities / hsp.align_length * 100 >= minimumIdentity):
                        annotation_status = check_annotation(f"{os.getcwd()}/bradyrhizobia-data", f"{alignment.hit_id}", f"{hsp.sbjct_start}", f"{hsp.sbjct_end}")
                        out_handle.write(f"{query.query}\t{alignment.hit_id}\t{hsp.identities / hsp.align_length * 100:.2f}\t{hsp.align_length / query.query_length * 100:.2f}\t{hsp.expect}\t{annotation_status}\n")
        result_handle.close()
    print(f"Filtering done. Check out {output_file} where the filtered results are saved")

def check_annotation(db_directory, hit_id, hit_start, hit_end):
    hit_id = hit_id.split("|")[1] #Example: NZ_CP029603.1
    assembly_id = get_assembly_identifier(hit_id) #Example: GCF_003183845.1
    locus_id = hit_id.split(".")[0] #Example: NZ_CP029603
    
    annotation_status = ""
    genbank_file_name = f"{assembly_id}_genomic.gbff"
    genbank_file_path = os.path.join(db_directory, genbank_file_name)
    if os.path.exists(genbank_file_path):
        for record in SeqIO.parse(genbank_file_path, "genbank"):
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

def get_assembly_identifier(hit_id):
    print(hit_id)
    handle = Entrez.efetch(id=hit_id, db='nucleotide', rettype='gb', retmode='text')
    efetch_record = SeqIO.read(handle, 'gb')
    db_cross_references = efetch_record.dbxrefs
    assembly = None
    biosample = None
    for reference in db_cross_references:
        if reference.startswith('Assembly:'):
            assembly = reference.split(':')[1]
        elif reference.startswith('BioSample:'):
            biosample = reference.split(':')[1]
        elif reference.startswith('BioProject:'):
            bioproject = reference.split(':')[1]
    
    if assembly:
        es_handle = Entrez.esearch(db='assembly', term=assembly)
        esearch_record = Entrez.read(es_handle)
    elif biosample:
        es_handle = Entrez.esearch(db='assembly', term=biosample)
        esearch_record = Entrez.read(es_handle)
    else:
        es_handle = Entrez.esearch(db='assembly', term=bioproject)
        esearch_record = Entrez.read(es_handle)
    
    
    for id in esearch_record['IdList']:
        es_handle = Entrez.esummary(db='assembly', id=id, report='full')
        es_record = Entrez.read(es_handle)
        document_summary = es_record["DocumentSummarySet"]["DocumentSummary"][0]
        repos = ["RefSeq", "GenBank"]
        for repo in repos:
            accession = document_summary.get(f"FtpPath_{repo}", False).split("/")[-1]
            if accession:
                return accession
        return None

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
    query = f"{os.getcwd()}/nodABC.fasta"
    db_name = f"{os.getcwd()}/bradyrhizobia-data/bradyrhizobia_blastdb"
    
    output_dir = "bradyrhizobia-blast"
    os.makedirs(output_dir, exist_ok=True)
    output = f"{os.getcwd()}/{output_dir}/tblastnResult_nodABC.xml"
    output_filtered = f"{os.getcwd()}/{output_dir}/tblastnResult_nodABC_Filtered"
    
    minimumIdentity = 30
    minimumCoverage = 50
    maximumEValue = 1e-5
    
    #blast_data(query, db_name, output)
    filter_blast_with_retry(output, output_filtered, minimumIdentity, minimumCoverage, maximumEValue)
