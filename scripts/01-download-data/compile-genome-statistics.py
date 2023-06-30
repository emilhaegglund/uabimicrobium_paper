"""
Script to merge the genome metadata downloaded from NCBI with the predefined groups
"""
import pandas as pd
import sys

accessions = pd.read_csv(sys.argv[1])
genome_metadata = pd.read_csv(sys.argv[2], sep="\t")

df = pd.merge(left=genome_metadata, right=accessions, on="Assembly Accession")
df.to_csv(sys.argv[3], sep="\t", index=False)
