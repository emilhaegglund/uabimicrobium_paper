# make protein kinase evaluation table of hmm proteins - for more realistic pkinase domain proteins.
import pandas as pd
import argparse

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument("--uab_hmm_file", required=True),
    parser.add_argument("--interproscan_file", required=True),
    parser.add_argument("--p_mapp_pvc", required=True), 
    parser.add_argument("--genome_stat", required=True),
    parser.add_argument("--out_file", required=True),
    parser.add_argument("--Analysis", required=True),
    parser.add_argument("--Domain", required=True),
    parser.add_argument("--uab_acc", required=True)
    return parser.parse_args()

def read_hmm_result(hmm_result_file):
	col_names = ['target name','accession_1', 'query name', 'accession_2', 'E-value_1', 'score_1', 'bias_1', 'E-value_2', 'score_2', 'bias_2', 'exp', 'reg', 'clu', 'ov', 'env', 'dom','rep', 'inc' ,'description of target']
	og_pkinas_df = pd.read_csv(hmm_result_file, names = col_names, sep='\t', comment='#')
	return og_pkinas_df

def open_interproscan_file(interproscan_file, analysis, selected_domain):
    """
    extract the Pfam domains of a interproscan file
    """
    interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession', 
                           'Signature description','Start location', 'Stop location', 'Score', 'Status', 
                           'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description', 
                           'optional: GO terms', 'optional: Pathway annotation']

    domain_df = pd.read_csv(interproscan_file, sep='\t', header=None, names = interproscan_header)
    domain_df = domain_df.fillna('Missing')
    domain_df = domain_df[domain_df['Analysis']==analysis]
    #print(domain_df['Signature accession'])
    pids = domain_df[domain_df['Signature accession']==str(selected_domain)]['protein_id'].values.tolist()
    #print(pids)
    domains_df = domain_df[domain_df['protein_id'].isin(pids)]
    
    domain_df = domain_df.sort_values(by=['protein_id','Start location'])
    
    return domain_df

def main(uab_hmm_file, interproscan_file, p_mapp_pvc, Analysis, domain, genome_stat, out_file, uab_acc):

    # open hmm file of uab and extract the protein id and evalue columns of the pkinase domains
    h_df = read_hmm_result(uab_hmm_file)
    h_df_pkinase_dom = h_df[['target name', 'E-value_1', 'E-value_2']]
    h_df_pkinase_dom.columns = ['protein_id', 'E-value_1', 'E-value_2']

    # open interprocasn file extract protein kinases and all their other domains
    i_df = open_interproscan_file(interproscan_file, Analysis, domain)

    protein_mapping = pd.read_csv(p_mapp_pvc, sep='\t')
    genome_stat = pd.read_csv(genome_stat, sep='\t')
    uab_pids = protein_mapping[protein_mapping['accession']==uab_acc]['protein_id'].drop_duplicates().values.tolist()
    print(uab_pids)
    print(i_df)
    # extract only the uab protein kinases and only the protein kinase domain
    #i_df_uab = i_df[i_df['protein_id'].isin(uab_pids)]['protein_id'].drop_duplicates().values.tolist()
    i_df_uab_df = i_df[i_df['protein_id'].isin(uab_pids)]
    print(i_df_uab_df)
    i_df_uab_df = i_df_uab_df[i_df_uab_df['Signature accession']==domain][['protein_id','Start location', 'Stop location', 'Score']]
    
    
    # merge the hmm info and interproscan info and save to file (also add the new column region length)
    #print('hmm')
    #print(h_df_pkinase_dom)
    #print('interpro')
    #print(i_df_uab_df)
    info_table = h_df_pkinase_dom.merge(i_df_uab_df, how='left', on = 'protein_id')
    
    info_table['region length'] = info_table['Stop location'] - info_table['Start location'] 
    info_table.columns = ['Protein id','E-value_1 hmm','E-value_2 hmm','Start location', 'Stop location', 'E-value interproscan', 'Region length']
    info_table.to_csv(out_file, sep='\t', index=False)
    return info_table


# input
args = parse_command_line()

uab_hmm_file = args.uab_hmm_file
interproscan_file = args.interproscan_file
p_mapp_pvc = args.p_mapp_pvc
genome_stat = args.genome_stat
out_file = args.out_file
Analysis = args.Analysis
Domain = args.Domain
uab_acc = args.uab_acc

info_table = main(uab_hmm_file, interproscan_file, p_mapp_pvc, Analysis, Domain, genome_stat, out_file, uab_acc)