import pandas as pd
import re
import os


def reshape_table(f):
    table = []
    temp = []

    for line in f:

        l = re.split(r'\s{2,}', line.rstrip())

        if l[0]=='ID':
            prot_id = l[1]

        if l[0] == 'FT':
            start = l[2] 
            stop =  l[3]

            if l[1] == 'SIGNAL':
                region_type = 'SIGNAL'
            if l[1] == 'TRANSMEM':
                region_type = 'TRANSMEM'
            if l[1] == 'DOMAIN':
                region_type = l[4]

            temp.append([prot_id, start, stop, region_type])

        if l[0] == '//':
            table+=temp
            temp=[]
    df = pd.DataFrame(table)
    df.columns = ['protein_id', 'start', 'stop', 'region_type']
    return df

    
def add_tm_region_type(phobius_df):
    new = [] 
    df_index = phobius_df[phobius_df['region_type'].isin(['SIGNAL','TRANSMEM'])].index
    phobius_df['region_type_detailed'] = phobius_df['region_type']
    for i in df_index:
        before = phobius_df.iloc[i-1]['region_type']
        after = phobius_df.iloc[i+1]['region_type']
        region = phobius_df.iloc[i]['region_type']
        
        if before == 'NON CYTOPLASMIC.':
            before = 'n'
        if before == 'CYTOPLASMIC.':
            before = 'c'
        if after == 'NON CYTOPLASMIC.':
            after = 'n'
        if after == 'CYTOPLASMIC.':
            after = 'c'
        if region == 'TRANSMEM':
            region = 'TM'
        
        if region != 'SIGNAL':
            region_type_added = before + '_' + region + '_' + after
        else:
            region_type_added = region 
            
        phobius_df.region_type_detailed.iloc[i] = region_type_added
    
    phobius_df = phobius_df[~phobius_df['region_type'].isin(['CYTOPLASMIC.', 'NON CYTOPLASMIC.', 'N-REGION.', 'H-REGION.', 'C-REGION.'])]
    return phobius_df

def main(phobus_file, df_info, phobius_reformat_out, tm_protein_table_out, tm_species_table_out):
    
    # phobius reformated table
    global phobius_df
    phobius_df = reshape_table(phobius_file)
    
    phobius_df_tm_type = add_tm_region_type(phobius_df)
    phobius_df_tm_type_out = phobius_df_tm_type[['protein_id', 'start', 'stop', 'region_type_detailed']]
    phobius_df_tm_type_out.columns = ['protein_id', 'start', 'stop', 'region_type']
    f = phobius_df_tm_type_out.to_csv(phobius_reformat_out, sep='\t', index=False)
    
    transmem_column = []
    signal_column = []
    cytoplasmic_column = []
    non_cytoplasmic_column = []
    protein_tm_table = []
    
    species_list = df_info['species'].values.tolist()
    
    for specie in species_list:
        pids = df_info[df_info['species']== specie]['pkinase protein ids'].values[0].split(';')
        annotation = df_info[df_info['species']== specie]['annotation'].values[0].split(';')

        transmem_counter = 0
        signal_counter = 0
        cytoplasmic_counter = 0
        non_cytoplasmic_counter = 0
        pid_amount = len(pids)
        counter = 0
        count_TM = 0
        
        for pid in pids:
            if pid in phobius_df['protein_id'].values.tolist():
                region_types = phobius_df[phobius_df['protein_id']==pid]['region_type'].values.tolist()
                
                if 'TRANSMEM' in region_types:
                    transmem_counter +=1
                        
                    a = annotation[counter]
                    count_TM = 0
                    if 'OS' in a:
                        a = a.split('OS')[0]
                    for r in region_types:
                        if 'TRANSMEM' == r:
                            count_TM += 1
                    protein_tm_table.append([specie, pid, a, 'TM', count_TM])
                
                if 'SIGNAL' in region_types:
                    signal_counter += 1
                
                if 'TRANSMEM' not in region_types: 
                        
                    a = annotation[counter]
                    if 'OS' in a:
                        a = a.split('OS')[0]
                    protein_tm_table.append([specie, pid, a, 'no TM', 0])
            else:
                protein_tm_table.append([specie, pid, annotation[counter], 'no TM', 0])

        
        transmem_column.append([specie, transmem_counter, round((int(transmem_counter)/int(pid_amount))*100)])
        signal_column.append([specie, signal_counter, round((int(signal_counter)/int(pid_amount))*100)])
    
    # save a table with per protein TM count
    protein_tm_table = pd.DataFrame(protein_tm_table)
    protein_tm_table.columns = ['specie','protein_id','annotation', 'TM_or_not', 'TM_count']
    protein_tm_table.to_csv(tm_protein_table_out, index=False, sep='\t')

    # create a table with the number of TM and signal proteins per species
    tm_species_table = df_info[['accession', 'species', 'group', 'nr of pkinases', 'all proteins']]
    
    transmem_df = pd.DataFrame(transmem_column)
    transmem_df.columns = ['species', 'transmembrane proteins', '% transmembrane proteins']
    tm_species_table = tm_species_table.merge(transmem_df, how='left', on='species')
    
    signal_df = pd.DataFrame(signal_column)
    signal_df.columns = ['species', 'signal proteins', '% signal proteins']
    tm_species_table = tm_species_table.merge(signal_df, how='left', on='species')
    # add species with zero tm regions
    acc_present = tm_species_table.accession.drop_duplicates().values.tolist()
    acc_abscent = tm_species_table[~tm_species_table['accession'].isin(acc_present)].drop_duplicates()
    acc_abscent['transmembrane proteins'] = 0
    acc_abscent['% transmembrane proteins'] = 0
    acc_abscent['signal proteins'] = 0
    acc_abscent['% signal proteins'] = 0 
    tm_species_table = pd.concat([tm_species_table, acc_abscent], ignore_index=True)
    s = tm_species_table.to_csv(tm_species_table_out, sep='\t', index=False)
    return tm_species_table, protein_tm_table

phobius_file = open(snakemake.input.phobius_file, "r") 
df_info = pd.read_csv(snakemake.input.protein_kinase_table, sep='\t')[['accession', 'species', 'group', 'nr of pkinases', 'all proteins', 'pkinase protein ids', 'annotation']] 
phobius_reformat_out = snakemake.output.phobius_reformat
tm_protein_table_out = snakemake.output.tm_protein_table
tm_species_table_out = snakemake.output.tm_species_table

tm_species_table, protein_tm_table = main(phobius_file, df_info, phobius_reformat_out, tm_protein_table_out, tm_species_table_out)  