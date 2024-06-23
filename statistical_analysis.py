import numpy as np
import os

from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
from pathlib import Path

def statistical_analysis(total_with_t6ss, total_without_t6ss, count_with_t6ss_1, count_without_t6ss_1):
    count_with_t6ss_0 = total_with_t6ss - count_with_t6ss_1
    count_without_t6ss_0 = total_without_t6ss - count_without_t6ss_1
    
    table = [
        [count_with_t6ss_1, count_with_t6ss_0],
        [count_without_t6ss_1, count_without_t6ss_0]
    ]

    # Chi-squared and Fisher's tests
    try:
        chi2, p_chi2, _, _ = chi2_contingency(table)
    except ValueError:
        p_chi2 = None

    try:
        oddsratio, p_fisher = fisher_exact(table)
    except ValueError:
        p_fisher = None

    # T-test and Mann-Whitney U test
    with_t6ss = [1] * count_with_t6ss_1 + [0] * count_with_t6ss_0
    without_t6ss = [1] * count_without_t6ss_1 + [0] * count_without_t6ss_0
    
    try:
        t_stat, p_ttest = ttest_ind(with_t6ss, without_t6ss)
    except ValueError:
        p_ttest = None
        
    try:
        u_stat, p_mannwhitney = mannwhitneyu(with_t6ss, without_t6ss, alternative='two-sided')
    except ValueError:
        p_mannwhitney = None

    return p_chi2, p_fisher, p_ttest, p_mannwhitney

def significance(p_value):
    if p_value:
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
    else:
        return "error"

if __name__ == "__main__":
    directory = os.path.join(os.getcwd(), "statistical-analysis")
    for file in Path(directory).rglob("*.tsv"):
        if file.stem.endswith("_output"):
            continue
        
        output_lines = []
        with open(f"{file}", "r") as adaptor_table:
            header = adaptor_table.readline().strip()
            output_lines.append("\t"+ header + "\tchi2\tfisher\tt-student\tmannwhitney")
            
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
                p_chi2, p_fisher, p_ttest, p_mannwhitney = statistical_analysis(total_with_t6ss,
                                                                                total_without_t6ss,
                                                                                nadaptors_t6ss,
                                                                                nadaptors_no_t6ss,)
                print(f"{adaptor}")
                print(f"- Chi-squared test p-value: {p_chi2}")
                print(f"- Fisher's exact test p-value: {p_fisher}")
                print(f"- T-test p-value: {p_ttest}")
                print(f"- Mann-Whitney U test p-value: {p_mannwhitney}\n")
                chi2_significance = significance(p_chi2)
                fisher_significance = significance(p_fisher)
                ttest_significance = significance(p_ttest)
                mannwhitney_significance = significance(p_mannwhitney)
                #output_lines.append(f"{adaptor}\t{nadaptors_t6ss}\t{nadaptors_no_t6ss}\t{chi2_significance}\t{fisher_significance}\t{ttest_significance}\t{mannwhitney_significance}")
                output_lines.append(f"{adaptor}\t{nadaptors_t6ss}\t{nadaptors_no_t6ss}\t{p_chi2}\t{p_fisher}\t{p_ttest}\t{p_mannwhitney}")

        output_file = f"{directory}/{file.stem}_output.tsv"
        with open(output_file, "w") as out_f:
            out_f.write("\n".join(output_lines))
                
