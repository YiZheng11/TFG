import os

def count1(filename, organism):
    organisms = {}
    organism_names = set()
    
    with open(filename, 'r') as file:
        file.readline()
        for line in file:
            columns = line.strip().split('\t')
            organism_def = columns[2]
            query_id = columns[0]
            
            if organism in organism_def:
                if organism_def not in organisms.keys():
                    organisms[organism_def] = [query_id]
                else:
                    if query_id not in organisms[organism_def]:
                        organisms[organism_def].append(query_id)

        for key in organisms.keys():
            if len(organisms[key]) >= 1:
                organism_names.add(key)
                    
    return len(organism_names), organism_names

def count3(filename, organism, queries1, queries2, queries3):
    organisms = {}
    organism_names = set()
    
    with open(filename, 'r') as file:
        file.readline()
        for line in file:
            columns = line.strip().split('\t')
            organism_def = columns[2]
            query_id = columns[0]
            
            if organism in organism_def:
                if organism_def not in organisms.keys():
                    organisms[organism_def] = [query_id]
                else:
                    if query_id not in organisms[organism_def]:
                        organisms[organism_def].append(query_id)
        
        for key in organisms.keys():
            found = False
            for query1 in queries1:
                for query2 in queries2:
                    for query3 in queries3:
                        if query1 in organisms[key]:
                            if query2 in organisms[key]:
                                if query3 in organisms[key]:
                                    organism_names.add(key)
                                    found = True
                                    break
                    if found:
                        break
                if found:
                    break
                    
    return len(organism_names), organism_names

def table_organisms(set1, set2):
    #unique_in_set1 = set1 - set2
    #unique_in_set2 = set2 - set1
    nodABC_nifH = set2.intersection(set1)

    #with open("both_rhizobium_nod_nif", 'w') as f1_unique:
    #    for organism in nodABC_nifH:
    #        f1_unique.write(organism + '\n')

    #with open("unique_nod", 'w') as f2_unique:
    #    for organism in unique_in_set1:
    #        f2_unique.write(organism + '\n')
    
    #with open("unique_nif", 'w') as f3_unique:
    #    for organism in unique_in_set2:
    #        f3_unique.write(organism + '\n')
    

    return len(nodABC_nifH)


if __name__ == "__main__":
    nodABC_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_nodABC_Filtered_noAnnotation"
    nodC_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_nodC_Filtered_noAnnotation"
    nifH_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_nifH_Filtered_noAnnotation"
    T6SSIII_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_T6SS_III_Filtered_noAnnotation"
    T6SSI_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_T6SS_I_Filtered_noAnnotation"
    T6SS_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_T6SS_Filtered_noAnnotation"
    adaptors_file = f"{os.getcwd()}/rizobia-blast/tblastnResult_PRK06147_Filtered_noAnnotation"
    
    queries_nodA = ["WP_057853260.1 NodA family N-acyltransferase [Bradyrhizobium valentinum]",
                    "WP_020923595.1 NodA family N-acyltransferase [Rhizobium etli]",
                    "WP_011084822.1 MULTISPECIES: NodA family N-acyltransferase [Bradyrhizobium]",
                    "WP_069091160.1 NodA family N-acyltransferase [Mesorhizobium sp. SEMIA 3007]",
                    "WP_015633488.1 MULTISPECIES: nodulation N-acyltransferase NodA [Sinorhizobium]",
                    "WP_018068349.1 MULTISPECIES: NodA family N-acyltransferase [Rhizobium]",
                    "WP_010967455.1 MULTISPECIES: nodulation N-acyltransferase NodA [Sinorhizobium]"]
    queries_nodB = ["WP_057853261.1 MULTISPECIES: chitooligosaccharide deacetylase NodB [Bradyrhizobium]",
                    "WP_020923594.1 chitooligosaccharide deacetylase NodB [Rhizobium etli]",
                    "WP_011084823.1 MULTISPECIES: chitooligosaccharide deacetylase NodB [Bradyrhizobium]",
                    "WP_069091164.1 MULTISPECIES: chitooligosaccharide deacetylase NodB [unclassified Mesorhizobium]",
                    "WP_010875356.1 MULTISPECIES: chitooligosaccharide deacetylase NodB [Sinorhizobium]",
                    "WP_018517097.1 chitooligosaccharide deacetylase NodB [Rhizobium leguminosarum]",
                    "WP_003532851.1 MULTISPECIES: chitooligosaccharide deacetylase NodB [Sinorhizobium]"]
    queries_nodC = ["WP_057853262.1 chitooligosaccharide synthase NodC [Bradyrhizobium valentinum]",
                   "WP_020923593.1 chitooligosaccharide synthase NodC [Rhizobium etli]",
                   "WP_011084824.1 MULTISPECIES: chitooligosaccharide synthase NodC [Bradyrhizobium]",
                   "WP_069091161.1 chitooligosaccharide synthase NodC [Mesorhizobium sp. SEMIA 3007]",
                   "WP_014858070.1 MULTISPECIES: chitooligosaccharide synthase NodC [Sinorhizobium]",
                   "WP_018517098.1 chitooligosaccharide synthase NodC [Rhizobium leguminosarum]",
                   "WP_010967454.1 chitooligosaccharide synthase NodC [Sinorhizobium meliloti]"]
    
    
    genera = ["Rhizobium", "Bradyrhizobium", "Mesorhizobium", "Sinorhizobium"]
    for genus in genera:
        #num1, lista1 = count3(nodABC_file, genus, queries_nodA, queries_nodB, queries_nodC)
        num2, lista2 = count1(adaptors_file, genus)
        print(genus, num2)
        #print(num2)
        #print("nodABC y nifH:", table_organisms(lista1, lista2))

