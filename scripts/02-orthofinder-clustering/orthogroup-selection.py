#!/user/bin/env/ python3

"""
Single copy gene clusters with at least 95% of the taxa present.
"""

import pandas as pd
import sys
import matplotlib.pyplot as plt


def extract_single_copy_OG(gene_count_table, OG_threshold):
    """
    extract single copy Orthogroups with chosen threshold of percentage of genomes present in the OG:s
    """
    count_genomes = len(gene_count_table.columns) - 2

    min_thresh = float(count_genomes) * float(OG_threshold)
    out_list = []
    og_all = []

    for i in gene_count_table.index:
        total = float(gene_count_table["Total"][i])
        OG = gene_count_table["Orthogroup"][i]
        og_all.append([total, OG])
        # check that the at least 95% of the genomes are present and that each genome only has one or less gene in the Orthogroup
        if (
            total >= min_thresh
            and all(gene_count < 2 for gene_count in gene_count_table.iloc[i, 1:-2])
            == True
        ):
            out_list.append([total, OG])
    return out_list, og_all


def save_to_file(filtered_OG, output_file_path):
    """
    save result to file and return the output file
    """
    with open(output_file_path, "w") as f:
        for line in filtered_OG:
            f.write("\t".join(map(str, line)) + "\n")


### Run the funcitons

# input
input_file_path = sys.argv[
    1
]  # path to Orthogroups.GeneCount.tsv produced by orthofinder
output_file_path = sys.argv[
    2
]  # the path to the output folder to store the selected Orthogroups file produced by this script

# 1. Open gene count file as a pandas dataframe
df = pd.read_csv(input_file_path, sep="\t")

# 2. Exctract the Orthogroup ids containing 95% or more of the genomes. (No more than one gene per genome is aloud)
filtered_OG, og_all = extract_single_copy_OG(df, 0.95)

# 3. Save the filtered total gene count and the Orthogroup id to file
save_to_file(filtered_OG, output_file_path)
