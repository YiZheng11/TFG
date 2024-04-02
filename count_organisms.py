import os

def count(filename, num_queries):
    organisms = {}
    organism_names = set()
    count = 0
    
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            organism_id = columns[1]
            query_id = columns[0]
            if organism_id not in organisms.keys():
                    organisms[organism_id] = [query_id]
            else:
                if query_id not in organisms[organism_id]:
                    organisms[organism_id].append(query_id)
        
        for key in organisms.keys():
            if len(organisms[key]) == num_queries:
                organism_names.add(key)
                count += 1
    return count, organism_names

def get_unique_organisms(set1, set2, output_name1, output_name2):
    unique_in_set1 = set1 - set2
    unique_in_set2= set2 - set1

    with open(output_name1, 'w') as f1_unique:
        for organism in unique_in_set1:
            f1_unique.write(organism + '\n')

    with open(output_name2, 'w') as f2_unique:
        for organism in unique_in_set2:
            f2_unique.write(organism + '\n')

    return None


if __name__ == "__main__":
    filename1 = f"{os.getcwd()}/bradyrhizobia-blast/tblastnResult_nodABC_Filtered_noAnnotation"
    filename2 = f"{os.getcwd()}/bradyrhizobia-blast/tblastnResult_nodC_Filtered_noAnnotation"
    output_name1 = f"{os.getcwd()}/bradyrhizobia-blast/nodABC_unique"
    output_name2 = f"{os.getcwd()}/bradyrhizobia-blast/nodC_unique"
    
    num_nodABC, list_nodABC = count(filename1, 3)
    #num_queries = 3 when using the file for nodABC
    num_nodC, list_nodC = count(filename2, 1)
    #num_queries = 1 when using the file for nodC
    print("Total number of organisms\n-nodABC:", num_nodABC, "\n-nodC:", num_nodC)
    
    #get_unique_organisms(list_nodABC, list_nodC, output_name1, output_name2)
