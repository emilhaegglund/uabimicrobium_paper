# from the proteome _ files create a file with GCF_acc and protein id column

from Bio import SeqIO
import pandas as pd
import os

protein_mapping = []
p_count_list = []

for fasta_file in os.listdir(snakemake.params):
    if fasta_file.split('.faa')[-1]!='.gz' and fasta_file.split('.faa')[0] in ['GCA_000002945.2','GCA_000146045.2','GCA_000002985.3','GCA_000005845.2']:
            group = ''
            species = ''
            
            if fasta_file.split('.faa')[0] == 'GCA_000002985.3':
                group = 'Nematoda'
                species = 'Caenorhabditis elegans'
                
            if fasta_file.split('.faa')[0] == 'GCA_000002945.2':
                group = 'Ascomycota'
                species = 'Schizosaccharomyces pombe'
                
            if fasta_file.split('.faa')[0] == 'GCA_000146045.2':
                group = 'Ascomycota'
                species = 'Saccharomyces cerevisiae'
                
            if fasta_file.split('.faa')[0] == 'GCA_000005845.2':
                group = 'Proteobacteria'
                species = 'Escherichia coli'
            
            p_count = 0
            
            
            for record in SeqIO.parse(snakemake.params + fasta_file ,'fasta'):
                protein_mapping.append([fasta_file.split('.faa')[0], record.id, group, species])
                p_count+=1
                
            p_count_list.append([fasta_file.split('.faa')[0], int(p_count)])
            
            
                

protein_mapping_df = pd.DataFrame(protein_mapping)
protein_mapping_df.columns = ['accession','protein_id', 'group', 'organism name']

for row in p_count_list:
    acc = row[0]
    p_count = row[1]
    protein_mapping_df.loc[protein_mapping_df[protein_mapping_df['accession']==acc].index, 'protein count'] = p_count

protein_mapping_df.to_csv(snakemake.output, sep='\t', index=False)
