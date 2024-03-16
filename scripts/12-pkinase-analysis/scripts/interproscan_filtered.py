import pandas as pd

def read_interpro(interpro_path):
    col_names = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession', 'Signature description', 'Start location', 'Stop location', 'Score', 'Status',
        'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
        'optional: GO terms', 'optional: Pathway annotation']
    interproscan_df = pd.read_csv(interpro_path, sep='\t', header=None, names=col_names)
    return interproscan_df


def main(interproscan_raw, pkinase_table, outfile):
    i_df = read_interpro(str(interproscan_raw))[['protein_id','Sequence length','Analysis', 'Signature accession', 'Signature description', 'Start location', 'Stop location', 'Score','Date of run']]
    
    pkinase_active_df = pd.read_csv(str(pkinase_table), sep='\t')
    pkinase_active_pid = pkinase_active_df['protein_id'].drop_duplicates().values.tolist()
    
    i_filtered_df = i_df[i_df['protein_id'].isin(pkinase_active_pid)]
    
    # replace or remove the pkinases that has suffered from insertions and deletions
    pkinases_to_change = pkinase_active_df[['protein_id','protein_from','protein_to']]#######OBS![pkinase_active_df['domain_completness']=='full merged'][['protein_id','protein_from','protein_to']]
    pkinases_to_change.columns = ['protein_id','Start location', 'Stop location']
    pkinases_to_change['Sequence length'] = ''
    pkinases_to_change['Analysis'] = 'Pfam'
    pkinases_to_change['Signature accession'] = 'PF00069'
    pkinases_to_change['Signature description'] = 'Protein kinase'
    pkinases_to_change['Score'] = ''
    pkinases_to_change['Date of run'] = ''
    
    i_filtered_df = i_filtered_df[i_filtered_df['Signature accession']!='PF00069']
    print(i_filtered_df)
    print(pkinases_to_change)
    i_filtered_df = pd.concat([i_filtered_df, pkinases_to_change], ignore_index=True)
    
    f = i_filtered_df.to_csv(str(outfile), index=False, sep='\t')
    
    return f

### run the script:

#input
interproscan_raw = str(snakemake.input.interproscan_raw)
pkinase_table = str(snakemake.input.pkinase_table)
outfile = str(snakemake.output.filtered) 

f = main(interproscan_raw, pkinase_table, outfile)