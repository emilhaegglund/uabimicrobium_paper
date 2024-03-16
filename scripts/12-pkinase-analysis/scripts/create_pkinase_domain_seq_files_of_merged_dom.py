from Bio import SeqIO
import pandas as pd
import argparse

# made emils script work for more than one species input file with functions and different input files

def create_domain_dict(df):
    
    domain_dict = {}

    for line_index in df.index:
        line_index = int(line_index)

        dom_id = df.iloc[line_index]['protein_id']
        start = df.iloc[line_index]['protein_from']
        stop = df.iloc[line_index]['protein_to']

        if dom_id not in domain_dict.keys():
            domain_dict[dom_id] = {"start": [start], "stop": [stop]}
        else:
            domain_dict[dom_id]["start"].append(start)
            domain_dict[dom_id]["stop"].append(stop)

    return domain_dict
            
def create_domain_seq_fasta(species_fasta_file, species_domain_out_fasta_file, domain_dict):
    with open(species_domain_out_fasta_file, "w") as f_out:
        for record in SeqIO.parse(species_fasta_file, "fasta"):
            if record.id in domain_dict.keys():
                for i in range(len(domain_dict[record.id]["start"])):
                    f_out.write('>' + record.id + '_pkinase_domain_' + str(i+1))
                    f_out.write('\n')
                    f_out.write(str(record.seq[domain_dict[record.id]["start"][i] : domain_dict[record.id]["stop"][i]]))
                    f_out.write('\n')
    f_out.close()
    return f_out

def create_group_fasta(species_domain_out_fasta_file, group_domain_out_fasta_file_prefix, domain_info_df):

    groups = domain_info_df['group'].drop_duplicates().values.tolist()
    
    for group in groups:
        pids = domain_info_df[domain_info_df['group']==group]['protein_id'].values.tolist()      
        f_group_out = group_domain_out_fasta_file_prefix + '_' + str(group) + '.fa'
        f_group_fasta = []
        
        for record in SeqIO.parse(species_domain_out_fasta_file, "fasta"):
            pid = record.id.split('_')[0]
            if pid in pids:
                f_group_fasta.append(record)
        #print(f_group_fasta)
        f_group = SeqIO.write(f_group_fasta, f_group_out, "fasta")     
        
    return f_group

def main(domain_info_file, species_fasta_file, species_domain_out_fasta_file, group_domain_out_fasta_file_prefix):
    domain_info_df = pd.read_csv(domain_info_file, sep='\t')
    domain_dict = create_domain_dict(domain_info_df)
    f = create_domain_seq_fasta(species_fasta_file, species_domain_out_fasta_file, domain_dict)
    f_group = create_group_fasta(species_domain_out_fasta_file, group_domain_out_fasta_file_prefix, domain_info_df)
    return f

# input:

domain_info_file = snakemake.input.domain_info_file
species_fasta_file = snakemake.input.species_fasta_file
species_domain_out_fasta_file = snakemake.output.species_domain_out_fasta_file
group_domain_out_fasta_file_prefix = snakemake.params.group_domain_out_fasta_file_prefix

# run the script:
f = main(domain_info_file, species_fasta_file, species_domain_out_fasta_file, group_domain_out_fasta_file_prefix)