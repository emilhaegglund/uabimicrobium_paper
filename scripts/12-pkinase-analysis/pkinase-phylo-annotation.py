"""
Script to create the annotation file for the
phylogeny of the Pkinase domain in the PVC.
"""
import os
import sys
import pandas as pd
from Bio import SeqIO

# Read domain architecture tables
df = pd.read_csv(sys.argv[1], sep="\t")
df["protein ids"] = df["protein ids"].str.split(";")
df = df.explode("protein ids", ignore_index=True)
df.rename(
    columns={
        "protein ids": "taxa",
        "Architecture": "architecture",
        "Domain annotations": "domain_annotation",
        "Protein annotation": "protein_annotation",
    },
    inplace=True,
)

# Read proteomes
proteome = {"taxa": [], "accession": []}
for f in os.listdir(sys.argv[2]):
    f_path = os.path.join(sys.argv[2], f)
    accession = os.path.splitext(f)[0]
    for record in SeqIO.parse(f_path, "fasta"):
        proteome["taxa"].append(record.id)
        proteome["accession"].append(accession)

proteome_df = pd.DataFrame.from_dict(proteome)

df = pd.merge(left=df, right=proteome_df, on="taxa", how="left")

# Read genome information
genome_df = pd.read_csv(sys.argv[3], sep="\t")

df = pd.merge(
    left=df, right=genome_df, left_on="accession", right_on="Assembly Accession"
)

df["taxa_new"] = df["taxa"] + "_1"

# Read sequences in phylogeny to find sequences with two
# pkinases
two_pkinase_copies = []
for record in SeqIO.parse(sys.argv[4], "fasta"):
    if record.id.split("_")[-1] == "2":
        two_pkinase_copies.append(record.id.split("_")[0])

# Test for absence of annotations
all_records = []
for record in SeqIO.parse(sys.argv[4], "fasta"):
    all_records.append(record.id.split("_")[0])
all_records = set(all_records)
df_records = set(df["taxa"].tolist())
print(all_records - df_records)
# Create a new df with sequences that have two copies and change
# the taxa_new col to end with _2
two_copies_df = df[df["taxa"].isin(two_pkinase_copies)]
two_copies_df["taxa_new"] = two_copies_df["taxa"] + "_2"
#print(two_copies_df)

# Concatenate the two data frames
df = pd.concat([df, two_copies_df])
df["taxa"] = df["taxa_new"]
df["arch-annotation"] = (
    df["taxa"] + " : " + df["architecture"] + " : " + df["Organism Name"] + " : " + df["Group_y"] + " : " + df["domain_annotation"]
)
df = df[["taxa", "architecture", "domain_annotation", "arch-annotation", "accession", "Group_y"]]
man_data = {"taxa":['OQD45576.1_1', 'GAB61147.1_1', 'GAN33544.1_1', 'SOH04169.1_1', 'AKC82788.1_1', 'QBG46553.1_1'], "architecture":[], "domain_annotation":[], "arch-annotation":[], "accession":[]
df.to_csv(sys.argv[5], sep="\t", index=False)
