from Bio import SeqIO
import os
import pandas as pd


# create a merged proteome file of all the proteomes

p_merged = snakemake.output
sequences = []

for proteome in os.listdir('../../results/01_download_data/genbank_proteomes/'):
  proteome_file = SeqIO.parse(open('../../results/01_download_data/genbank_proteomes/' + str(proteome), 'r'), 'fasta')

  for protein in proteome_file:
      sequences.append(protein)

SeqIO.write(sequences, str(p_merged), 'fasta')
