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

    # Chi-squared and Fisher's exact tests
    try:
        chi2, p_chi2, _, _ = chi2_contingency(table)
    except ValueError as e:
        p_chi2 = f"{e}"

    try:
        oddsratio, p_fisher = fisher_exact(table)
    except ValueError as e:
        p_fisher = f"{e}"

    # T-test and Mann-Whitney U test
    with_t6ss = [1] * count_with_t6ss_1 + [0] * count_with_t6ss_0
    without_t6ss = [1] * count_without_t6ss_1 + [0] * count_without_t6ss_0
    
    t_stat, p_ttest = ttest_ind(with_t6ss, without_t6ss)
    u_stat, p_mannwhitney = mannwhitneyu(with_t6ss, without_t6ss, alternative='two-sided')

    return p_chi2, p_fisher, p_ttest, p_mannwhitney

if __name__ == "__main__":
    directory = os.path.join(os.getcwd(), "adaptors_table")
    for file in Path(directory).rglob("*.tsv"):
        with open(f"{file}", "r") as adaptor_table:
            header = adaptor_table.readline()
            
            total = adaptor_table.readline().split("\t")
            total_with_t6ss = int(total[1])
            total_without_t6ss = int(total[2])
            
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