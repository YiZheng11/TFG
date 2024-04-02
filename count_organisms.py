import os

def count(filename):
    organisms = set()
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            organism_id = columns[1]
            query_id = columns[0]
            if query_id.count(organism_id) >= 3:
                    organisms.add(organism_id)
    return len(organisms)

def get_unique_organisms(filename1, filename2, output_name1, output_name2):
    organisms_file1 = set()
    organisms_file2 = set()

    with open(filename1, 'r') as file1:
        for line in file1:
            columns = line.strip().split('\t')
            organisms_file1.add(columns[1])

    with open(filename2, 'r') as file2:
        for line in file2:
            columns = line.strip().split('\t')
            organisms_file2.add(columns[1])

    file1_unique = organisms_file1 - organisms_file2
    file2_unique= organisms_file2 - organisms_file1

    with open(output_name1, 'w') as f1_unique:
        for organism in file1_unique:
            f1_unique.write(organism + '\n')

    with open(output_name2, 'w') as f2_unique:
        for organism in file2_unique:
            f2_unique.write(organism + '\n')

    return None


if __name__ == "__main__":
    filename1 = f"{os.getcwd()}/bradyrhizobia-blast/tblastnResult_nodABC_Filtered_noAnnotation"
    filename2 = f"{os.getcwd()}/bradyrhizobia-blast/tblastnResult_nodC_Filtered_noAnnotation"
    output_name1 = f"{os.getcwd()}/bradyrhizobia-blast/nodABC_unique"
    output_name2 = f"{os.getcwd()}/bradyrhizobia-blast/nodC_unique"
    
    num_organisms = count(filename2)
    print("Total number of organisms:", num_organisms)
    
    #get_unique_organisms(filename1, filename2, output_name1, output_name2)
