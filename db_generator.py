import os
import urllib
import gzip
import pdb
import subprocess

from pathlib import Path

from Bio import Entrez
from Bio import SeqIO

Entrez.email = "yi.zheng@alumnos.upm.es" # your email please

def get_index(number):
    "Return a NCBI summary of the assembly with id <number>"
    handle = Entrez.esummary(db="assembly", id=number, report="full")
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

    for number in record["IdList"]:
        index = get_index(number)
        url = get_url(index)
        
        # Check if URL is None
        if url is None:
            print(f"No download URL found for assembly {number}. Skipping...")
            continue
        
        name = Path(url).name
        files = [f"{url}/{name}_genomic.fna.gz", f"{url}/{name}_genomic.gbff.gz"]
        for filename in files:
            print(f"Fetching {filename} from NCBI")
            with urllib.request.urlopen(filename) as response:
                with open(f"{target_dir}/{Path(filename).name}", "wb") as outfile:
                    outfile.write(response.read())

    return None

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
    
    if fastas:
        os.system(f"makeblastdb -in '{fastas}' -title {title} -out {directory}/{title} -parse_seqids -dbtype nucl")
    return None

if __name__ == "__main__":
    data_dir = f"{os.getcwd()}/data"
    os.makedirs("data", exist_ok=True)
    get_assemblies("bradyrhizobium", data_dir, retmax=20)
    unzip_all(data_dir)
    directory_to_database(data_dir, "bradyrhizobia_blastdb")
