# make table of the KDD sites, extract the protein kinases with at least one domain with KDD site
# from the alignment extract positions with 90% shared K or D
import pandas as pd
from Bio import AlignIO
import os
#filename = "../results/active_sites/domain_msa/PF00069_Saltatorellus.aln"


def create_aa_site_df(filename):

    KD = {}
    K_pos = []
    D_pos = []
    pos = 0
    protein_amount = 0
    format = "fasta"
    
    for alignments in AlignIO.parse(filename, format):

        for protein in alignments:
            protein_amount+=1
            pos = 0

            for aa in protein.seq:
                pos+=1

                if str(aa) == 'K':
                    K_pos.append([protein.id, pos])

                if str(aa) == 'D':
                    D_pos.append([protein.id, pos])

        K_count = pd.DataFrame(K_pos, columns = ['protein id','position'])['position'].value_counts()
        K_df = K_count.to_frame().reset_index()
        K_df.columns = ['K pos.', 'count']
        K_df['count'].astype(int)
        K_df['%K'] = (K_df['count'] / int(protein_amount))*100
        K_df = K_df[K_df['%K']>60]#.values.tolist()

        D_count = pd.DataFrame(D_pos, columns = ['protein id','position'])['position'].value_counts()
        D_df = D_count.to_frame().reset_index()
        D_df.columns = ['D pos.', 'count']
        D_df['count'].astype(int)
        D_df['%D'] = (D_df['count'] / int(protein_amount))*100
        D_df = D_df[D_df['%D']>60]#.values.tolist()
        
    return K_df, D_df
    
def main(aln_folder, out_K, out_D):
    K = []
    D = []
    for aln_file in os.listdir(aln_folder):
        K_df, D_df = create_aa_site_df(aln_folder + '/' + aln_file)
        
        K_df['msa'] = aln_file
        D_df['msa'] = aln_file
        
        D+=D_df.values.tolist()
        K+=K_df.values.tolist()
    
    K = pd.DataFrame(K, columns=['K pos.','count', '%K', 'msa'])
    D = pd.DataFrame(D, columns=['D pos.','count', '%D', 'msa'])
    
    K.to_csv(out_K, index=False, sep='\t')
    D.to_csv(out_D, index=False, sep='\t')
    return K, D

k,d = main(snakemake.params.msa, snakemake.output.K, snakemake.output.D)
#k, d = main('../results/active_sites/domain_msa/')
