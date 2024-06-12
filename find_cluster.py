import os
import subprocess

from pathlib import Path

def makedb(files, dbname):
    if files:
        subprocess.Popen(f"makedb {dbname} {files}", shell=True).wait()

def multigeneblast(query_file, database, output):
    "Run MultiGeneBlast..."
    
    subprocess.Popen(f"multigeneblast -in {query_file} -db {database} -out {output}", shell=True).wait()

if __name__ == "__main__":
    files = "*.gbk"
    dbname = "rizobia_mgbdb"
    #makedb(files, dbname)
    
    query_file = os.path.join(os.getcwd(), "T6SS_mgb.fasta")
    
    for i in range(18):
        n = i + 1
        db_path = f"/home/yi/Documentos/TFG_Yi/rizobia-T6SS-mgb/rizobia-T6SS-mgb-{n}/rizobia_mgbdb"
        output_dir = f"/home/yi/Documentos/TFG_Yi/rizobia_cluster/rizobia_cluster_{n}"
        os.makedirs(f"{output_dir}", exist_ok=True)
        multigeneblast(query_file, db_path, output_dir)
        print(f"Run {n} MultiGeneBlast finished")
