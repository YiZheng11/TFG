import os
import urllib
import gzip
import pdb
import time
from urllib.error import HTTPError 
import subprocess

from pathlib import Path

from Bio import Entrez
from Bio.SeqIO import parse


Entrez.email = "lewis.grozinger@upm.es" # your email please
# Entrez.api_key = "1c105008de567a0fdcc74eedb9584b2ec109"  # your API key

def get_index(number):
    "Return a NCBI summary of the assembly with id <number>"
    handle = Entrez.esummary(db="assembly", id=number, report="full")
    return Entrez.read(handle, validate=False) #To skip all tags that are not represented in the DTD

def get_url(index):
    "Returns the url of the available download directory, trying each repo in <repos>"
    document_summary = index["DocumentSummarySet"]["DocumentSummary"]
    if document_summary:
        repos = ["RefSeq", "GenBank", "Assembly_rpt"]
        for repo in repos:
            url = document_summary[0].get(f"FtpPath_{repo}", False)
            if url:
                if url.endswith('.txt'):
                    print(f"URL '{url}' does not point to a directory. Skipping...")
                else:
                    return url
    return None

def get_assembly(number, target_dir, tries=1):
    index = get_index(number)
    url = get_url(index)
        
    if url is not None:
        while tries > 0:
            try: 
                name = Path(url).name
                file_url = f"{url}/{name}_genomic.gbff.gz"
                print(f"Fetching {file_url} from NCBI")
                output_fn = f"{target_dir}/{Path(file_url).name}"
                with urllib.request.urlopen(file_url) as response:
                    with open(output_fn, "wb") as outfile:
                        outfile.write(response.read())
                return output_fn
            except:
                print(f"Attempt for {file_url} failed, {tries} attempts remaining")
                time.sleep(3.33)
    else:
        print(f"No download URL found for assembly {number}. Skipping...")
        
    return None
    

def get_assemblies(term, target_dir, retmax=0):
    "Download to <target_dir> all fasta and genbank for assemblies matching <term>"
    handle = Entrez.esearch(db="assembly", term=term, retmax=retmax)
    record = Entrez.read(handle)
    N = len(record['IdList'])
    print(f"Found {N} records matching {term}")
    
    for (i, number) in enumerate(record["IdList"]):
        print(f"{i}/{N}")
        newfile = get_assembly(number, target_dir, tries=10)
        if newfile is not None:
            unzip(newfile)
        
    return None

def unzip(path_or_string):
    path = Path(path_or_string)
    print(f"Gunzipping {path}")

    with gzip.open(path, "rb") as zipped:
        with open(path.parent / path.stem, "wb") as unzipped:
            unzipped.write(zipped.read())
        
    path.unlink()

def unzip_all(base_dir):
    "Use gunzip to decompress all gzipped files under <base_dir>"
    base_dir = Path(base_dir)
    for fn in base_dir.rglob("*.gz"):
        unzip(base_dir / fn)
        
def directory_to_database(directory, title):
    "Recursively look for all fasta files in directory add build a blast db called <title> from them"
    directory = Path(directory)
    processed_ids = {}
    failed_ids = set()
    
    command = ["makeblastdb"]
    command += ["-in", "-"]
    command += ["-title", title]
    command += ["-out", f"{directory}/{title}"]
    command += ["-parse_seqids"]
    command += ["-dbtype", "nucl"]
    
    with subprocess.Popen(command, stdin=subprocess.PIPE) as process:
        for path in directory.rglob("*.gbff"):
            for record in parse(path, "genbank"):
                print(f"Adding {record.id} to database {title}")
                try:
                    if record.id not in processed_ids:
                        process.stdin.write(record.format("fasta").encode())
                        processed_ids[record.id] = 1
                    else:
                        processed_ids[record.id] += 1
                except:
                    print(f"Failed adding {record.id} to database {title}")
                    failed_ids.add(record.id)

    duplicates = "\n".join(filter(lambda x: processed_ids[x] > 1, processed_ids.keys()))
    print("Duplicates:")
    print(duplicates)
    print(f"Failed:")
    print(failed_ids)
    return None

def generate(term):
    data_dir = f"{os.getcwd()}/{term}-data"
    os.makedirs(f"{data_dir}", exist_ok=True)
    # get_assemblies(term, data_dir, retmax=2)
    directory_to_database(data_dir, f"{term}_blastdb")

if __name__ == "__main__":
    generate("Sinorhizobium")
