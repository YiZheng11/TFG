import numpy as np
import os

from scipy.stats import binomtest
from pathlib import Path

def binomial_proportion_test(total_with_t6ss, total_without_t6ss, count_with_t6ss_1, count_without_t6ss_1):
    total = total_with_t6ss + total_without_t6ss
    
    # Combined proportions
    combined_prop = (count_with_t6ss_1 + count_without_t6ss_1) / (total_with_t6ss + total_without_t6ss) if total > 0 else 0
    
    result = binomtest(count_with_t6ss_1, n=total_with_t6ss, p=combined_prop, alternative='two-sided')
    
    return result.pvalue

def significance(p_value):
    if p_value < 0.0001:
        return '****'
    elif p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''

if __name__ == "__main__":
    directory = os.path.join(os.getcwd(), "statistical-analysis")
    for file in Path(directory).rglob("*.tsv"):
        if file.stem.endswith("_output"):
            continue
        
        output_lines = []
        with open(f"{file}", "r") as adaptor_table:
            header = adaptor_table.readline().strip()
            output_lines.append("\t"+ header + "\tbinomial")
            
            total = adaptor_table.readline().strip().split("\t")
            total_with_t6ss = int(total[1])
            total_without_t6ss = int(total[2])
            output_lines.append("\t".join(total))
            
            print(file, "\n", total_with_t6ss, total_without_t6ss)
            for line in adaptor_table:
                adaptor_info = line.split("\t")
                adaptor = adaptor_info[0]
                nadaptors_t6ss = int(adaptor_info[1])
                nadaptors_no_t6ss = int(adaptor_info[2])
                
                p_binom = binomial_proportion_test(total_with_t6ss, total_without_t6ss, nadaptors_t6ss, nadaptors_no_t6ss)
                print(f"{adaptor}")
                print(f"- Binomial proportion test p-value: {p_binom}\n")
                binom_significance = significance(p_binom)
                output_lines.append(f"{adaptor}\t{nadaptors_t6ss}\t{nadaptors_no_t6ss}\t{binom_significance}")
                
        output_file = f"{directory}/{file.stem}_binom_output.tsv"
        with open(output_file, "w") as out_f:
            out_f.write("\n".join(output_lines))
                
