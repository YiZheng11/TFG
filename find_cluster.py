import os
import subprocess

def multigeneblast(query_file, database, start, end, genes):
    "Run MultiGeneBlast..."
    
    command = ["multigeneblast"]
    command += ["-in", query_file]
    command += ["-db", database]
    command += ["-from", start]
    command += ["-to", end]
    command += ["-genes", genes]

    subprocess.run(command)

if __name__ == "__main__":
    query_file = os.path.join(os.getcwd(), "GCF_000442435.1.gbff")
    database = os.path.join(os.getcwd(), "rizobia-data", "rizobia_blastdb")
    
    genes_accession = """WP_020919499.1, WP_020919500.1, WP_020919501.1, 
    WP_020919503.1, WP_020919504.1, WP_020919505.1, WP_020919507.1, 
    WP_041679311.1, WP_020919510.1, WP_041679313.1, WP_020919514.1"""
    
    start = "491148"
    end = "515517"

    multigeneblast(query_file, database, start, end, genes_accession)