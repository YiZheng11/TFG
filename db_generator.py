import os
import urllib
import gzip
import pdb
import time
from urllib.error import HTTPError 
#import subprocess

from pathlib import Path

from Bio import Entrez

Entrez.email = "YOUR_EMAIL" # your email please
Entrez.api_key = "YOUR_API_KEY"  # your API key

def get_index(number):
    "Return a NCBI summary of the assembly with id <number>"
    handle = Entrez.esummary(db="assembly", id=number, report="full")
    time.sleep(0.15)
    return Entrez.read(handle)

def get_url(index):
    "Returns the url of the available download directory, trying each repo in <repos>"
    document_summary = index["DocumentSummarySet"]["DocumentSummary"][0]
    repos = ["RefSeq", "GenBank", "Assembly_rpt"]

    for repo in repos:
        url = document_summary.get(f"FtpPath_{repo}", False)
        if url:
            return url

    return None

def get_assemblies(term, target_dir, retmax=0):
    "Download to <target_dir> all fasta and genbank for assemblies matching <term>"
    handle = Entrez.esearch(db="assembly", term=term, retmax=retmax)
    record = Entrez.read(handle)
    print(f"Found {len(record['IdList'])} records matching {term}")
    time.sleep(0.15)

    for number in record["IdList"]:
        index = get_index(number)
        url = get_url(index)
        
        # Check if URL is None
        if url is None:
            print(f"No download URL found for assembly {number}. Skipping...")
            continue
        
        name = Path(url).name
        files = [f"{url}/{name}_genomic.fna.gz", f"{url}/{name}_genomic.gbff.gz"]
        for file_url in files: 
            filename = Path(file_url).name
            filepath = Path(target_dir) / filename

            if not filepath.exists():
                print(f"Downloading {filename} from NCBI")
                with urllib.request.urlopen(url) as response:
                    with open(filepath, "wb") as outfile:
                        outfile.write(response.read())
            else: # In case some of the assemblies have already been downloaded
                print(f"Skipping {filename}, already exists")

    return None

def get_assemblies_with_retry(term, target_dir, retmax=0, max_retries=20):
    "Download to <target_dir> all fasta and genbank for assemblies matching <term>, with retry"
    retries = 0
    success = False

    while not success and retries < max_retries:
        try:
            get_assemblies(term, target_dir, retmax)
            success = True
        except HTTPError as e:
            if e.code == 500:  # Internal Server Error
                print("Internal Server Error, retrying again in 60 seconds...")
                time.sleep(60)
                retries += 1
            else:
                raise  # If it is not an Internal Server Error, it shows the exception/error captured

    if not success:
        print(f"Failed to download the assemblies for {term} after {max_retries} retries.")

def unzip_all(base_dir):
    "Use gunzip to decompress all gzipped files under <base_dir>"
    base_dir = Path(base_dir)
    for fn in base_dir.rglob("*.gz"):
        fn = Path(fn)
        print(f"Gunzipping {fn}")
        with gzip.open(fn, "rb") as zipped:
            with open(base_dir / fn.stem, "wb") as unzipped:
                unzipped.write(zipped.read())
        
        # Delete the .gz file after successfully unzipping
        fn.unlink()

def directory_to_database(directory, title):
    "Recursively look for all fasta files in directory add build a blast db called <title> from them"
    directory = Path(directory)
    fastas = " ".join(map(str, directory.rglob("*.fna")))
    
    #if fastas:
    #    command = f"makeblastdb -in '{fastas}' -title {title} -out {directory}/{title} -parse_seqids -dbtype nucl"
    #    if command is not None:
    #        subprocess.Popen(command, shell=True).wait()
    
    if fastas:
        os.system(f"makeblastdb -in '{fastas}' -title {title} -out {directory}/{title} -parse_seqids -dbtype nucl")
    return None

if __name__ == "__main__":
    data_dir = f"{os.getcwd()}/data"
    os.makedirs("data", exist_ok=True)
    
    organisms = ['Rhizobium', 'Bradyrhizobium', 'Mesorhizobium', 
                 'Sinorhizobium', 'Neorhizobium', 'Georhizobium', 
                 'Pararhizobium', 'Pseudorhizobium']
    
    for organism in organisms:
        get_assemblies_with_retry(organism, data_dir, retmax=100000)
    
    unzip_all(data_dir)
    directory_to_database(data_dir, "rizobia_blastdb")
