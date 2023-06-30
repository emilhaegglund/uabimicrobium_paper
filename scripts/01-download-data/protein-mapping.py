"""
Creates a mapping between protein identifier and species.
"""
from Bio import SeqIO
import sys
import os
import pandas as pd

mapping_dict = {}
for faa_path in sys.argv[1:-1]:
    # Extract the accession from the file name
    accession = os.path.split(faa_path)[1]
    accession = os.path.splitext(accession)[0]
    for record in SeqIO.parse(faa_path, 'fasta'):
        mapping_dict[record.id] = [accession]

df = pd.DataFrame.from_dict(mapping_dict)
df = df.T
df.reset_index(inplace=True)
df.rename(columns={0:'accession', 'index':'protein_id'}, inplace=True)
df.to_csv(sys.argv[-1], sep='\t', index=False)
