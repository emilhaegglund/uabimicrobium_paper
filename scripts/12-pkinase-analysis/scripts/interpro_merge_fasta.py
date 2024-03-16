import pandas as pd
import os
import sys
from Bio import SeqIO
import argparse

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteome_path", required=True),
    parser.add_argument("--out_file", required=True),
    parser.add_argument("--hmm_file", required=True)
    return parser.parse_args()

def read_hmm_result(hmm_result_file):
    col_names = ['target name','accession_1', 'query name', 'accession_2', 'E-value_1', 'score_1', 'bias_1', 'E-value_2', 'score_2', 'bias_2', 'exp', 'reg', 'clu', 'ov', 'env', 'dom','rep', 'inc' ,'description of target']
    og_pkinas_df = pd.read_csv(hmm_result_file, names = col_names, sep='\t', comment='#')
    return og_pkinas_df

#def remove_bad_evalues(df):
    #df = df[df['E-value_1'] < 1e-20]
    #df = df[df['E-value_2'] < 1e-20]
    #return df

def extract_pid(df):
    pid = df['target name'].drop_duplicates().values.tolist()
    return pid

def merge_proteomes(out_file, proteom_file_path, pid_list):
    
    #check if the output file exists!
    if os.path.isfile(out_file):
        merged_fasta_file = open(out_file, 'a')
    else:
        merged_fasta_file = open(out_file, 'w')
    
    # open the proteome file and write the extracted proteins to the outputfile 
    in_file = SeqIO.parse(open(proteom_file_path, 'rt'),'fasta')
    
    for sequence in in_file:
        if sequence.id in pid_list:
            SeqIO.write([sequence], merged_fasta_file, 'fasta')
    
    merged_fasta_file.close()
    return merged_fasta_file

  
def main(hmm_file, out_file, proteom_file_path):

    pid_list = []

    df = read_hmm_result(hmm_file)
    pids = extract_pid(df)
    f = merge_proteomes(str(out_file), str(proteom_file_path), pids)
    return f

# input
args = parse_command_line()

out_file = args.out_file
proteom_file_path = args.proteome_path
hmm_file = args.hmm_file

# Run the functions  
f = main(hmm_file, out_file, proteom_file_path) 


