#!/usr/bin/env python
import pandas as pd
import os
import re
import sys
import argparse

def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument('--protein_families', required=True)
    parser.add_argument('--out', required=True)

    return parser.parse_args()

args = parse_command_line()

# Read orthogroup gene count file
df = pd.read_csv(args.protein_families, sep='\t')

# Make the orthogroup names to column
#df.rename(columns={'Unnamed: 0':'Family'}, inplace=True)
#df.drop(columns='Total', inplace=True)

# Remove underscores from column names
print(list(df.columns.values)[0])
new_names = {}
for column_name in list(df.columns.values)[1:]:  # Skip the Orthogroup column
    print(column_name)
    accession = os.path.split(column_name)[1]
    print(accession)
    accession = re.search("GC(?:A|F)_[0-9]*\.[0-9]|gobs|cjuql4|soil9|mblw", accession)
    accession = accession.group()  # Use only the accession number as identifier.
    new_column_name = accession.replace('_', '')
    new_names[column_name] = new_column_name
df.rename(columns=new_names, inplace=True)

# Save file
df.to_csv(args.out, sep='\t', index=False)
