# Extract the nucleotide sequnces of selected orthogroups for the dN analysis.

from Bio import SeqIO
import sys
import os
import pandas as pd

# input
input_aa_path = snakemake.input[0] #'../../results/07_align_sco/trimmed_alignments/'
input_nt_path = snakemake.input[1] #'../../results/01_download_data/genbank_genes/'
protein_mapping = snakemake.input[2] #'../../results/01_download_data/protein_mapping.tsv'

# output
output_path = snakemake.output[0] #'../../results/10_species_phylogeny_slow_evolving/orthogroup_nt_sequences/'

# extract the nt sequences and save to files

def extract_nt_seq(input_nt_path, input_aa_path, output_file, protein_mapping):

    # variables
    folder_name = ''
    first_seq = True
    p_mapp = pd.read_csv(protein_mapping, sep='\t', header = None)

    # open each concatenated msa amino acid fasta file, one for each cluster
    aa_seq = SeqIO.parse(open(str(input_aa_path)), 'fasta')

    # for each sequence in the protein msa file search for the gene id in the .cds.fna file
    for seq in aa_seq:
        file_name = p_mapp[p_mapp.iloc[:,0]==seq.id].iloc[:,1].to_numpy()[0] # Use the protein_mapping file to extract the genome name for the protein id
        nt_seq = SeqIO.parse(os.path.join(input_nt_path, file_name + '.cds.fna'), 'fasta')

        # If the protein_id match extract the corresponding gene id save the nucleotide sequences to a new nucleotide OG file


        for s in nt_seq:
            if str(s.id.split('|')[-1].split('_')[2]).rstrip() == str(seq.id).rstrip():
                    if first_seq == True:
                        f = open(output_file, 'w')
                        SeqIO.write([s], f, "fasta")
                        first_seq = False
                    else:
                        f = open(output_file, 'a')
                        SeqIO.write([s], f, "fasta")

    return f


# Run the function

extract_nt_seq(input_nt_path, input_aa_path, output_path, protein_mapping)
