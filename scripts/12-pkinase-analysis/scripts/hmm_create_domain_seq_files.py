from Bio import SeqIO
import pandas as pd
import argparse

# made emils script work for more than one species input file with functions

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmm_domain_info_file", required=True),
    parser.add_argument("--species_fasta_file", required=True),
    parser.add_argument("--species_domain_out_fasta_file", required=True)
    return parser.parse_args()    

def create_domain_dict(hmm_domain_info_file):
    domain_dict = {}
    
    with open(hmm_domain_info_file, "r") as f:
        for line in f:
            if not line[0] == "#":
                line = line.strip()
                line = line.split()
                start = int(line[17])
                stop = int(line[18])
                if line[0] not in domain_dict.keys():
                    domain_dict[line[0]] = {"start": [start], "stop": [stop]}
                else:
                    domain_dict[line[0]]["start"].append(start)
                    domain_dict[line[0]]["stop"].append(stop)
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

def main(hmm_domain_info_file, species_fasta_file, species_domain_out_fasta_file):
    domain_dict = create_domain_dict(hmm_domain_info_file)
    f = create_domain_seq_fasta(species_fasta_file, species_domain_out_fasta_file, domain_dict)
    return f

# input:
args = parse_command_line()
hmm_domain_info_file = args.hmm_domain_info_file
species_fasta_file = args.species_fasta_file
species_domain_out_fasta_file = args.species_domain_out_fasta_file

# run the script:
f = main(hmm_domain_info_file, species_fasta_file, species_domain_out_fasta_file)