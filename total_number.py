import os
from Bio import SeqIO

def number_organism(directory):
    rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium = 0, 0, 0, 0
    number_files = 0
    for file in os.listdir(directory):
        if file.endswith(".fna"):
            number_files += 1
            file = os.path.join(directory, file)
            for record in SeqIO.parse(file, "fasta"):
                if "Rhizobium" in record.description:
                    rhizobium += 1
                    break
                elif "Bradyrhizobium" in record.description:
                    bradyrhizobium += 1
                    break
                elif "Mesorhizobium" in record.description:
                    mesorhizobium += 1
                    break
                elif "Sinorhizobium" in record.description:
                    sinorhizobium += 1
                    break
                else:
                    print(file)
                    break
    return number_files, rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium

if __name__ == "__main__":
    directory = f"{os.getcwd()}/rizobia-data"
    number_files, rhizobium, bradyrhizobium, mesorhizobium, sinorhizobium = number_organism(directory)
    print("Número total:", number_files, "\nRhizobium:", rhizobium,
          "\Bradyrhizobium:", bradyrhizobium, "\Mesorhizobium:", mesorhizobium,
          "\nSinorhizobium:", sinorhizobium)