
import argparse
import pandas as pd
from collections import Counter 

def open_interproscan_file(interproscan_file, analysis):
    """
    extract the Pfam domains of a interproscan file
    """
    #interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession', 
    #                       'Signature description','Start location', 'Stop location', 'Score', 'Status', 
    #                       'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description', 
    #                       'optional: GO terms', 'optional: Pathway annotation']
    interproscan_header = ['protein_id','Sequence length','Analysis', 'Signature accession', 'Signature description', 'Start location', 'Stop location', 'Score','Date of run']
    domain_df_raw = pd.read_csv(interproscan_file, sep='\t', header=None, names = interproscan_header)
    domain_df = domain_df_raw.fillna('Missing')
    domain_df = domain_df[domain_df['Analysis']==analysis]
    domain_df = domain_df.sort_values(by=['protein_id','Start location'])

    return domain_df

def create_df_info(protein_mapping, genome_stat, genome_stat_p_mapp_eu):
    p_mapp_df = pd.read_csv(protein_mapping, sep='\t')
    genome_stat_df = pd.read_csv(genome_stat, sep='\t')[['organism_name', 'group', 'assembly_accession']]
    genome_stat_df.columns = ['organism_name', 'group', 'accession']
    genome_stat_p_mapp_eu_df = pd.read_csv(genome_stat_p_mapp_eu, sep='\t')[['organism name', 'group', 'accession', 'protein_id']]
    genome_stat_p_mapp_eu_df.columns = ['organism_name', 'group', 'accession', 'protein_id']
    
    df_info = genome_stat_df.merge(p_mapp_df, how='left', on='accession')
    df_info = pd.concat([df_info,genome_stat_p_mapp_eu_df])[['organism_name', 'accession', 'group', 'protein_id']]
    return df_info

def create_dom_count(df, df_info):
    df = df.merge(df_info, how='left', on='protein_id')
    
    df_dom = df.drop(columns = ['protein_id'])
    df_dom = pd.DataFrame(df_dom.value_counts()).reset_index()

    df_dom.columns = ['Signature accession merged', 'Signature description merged', 'accession', 'organism_name', 'group', '#of domains']

    df_prot = df
    df_prot = df_prot.drop_duplicates() 
    df_prot = df_prot.drop(columns = ['protein_id'])
    df_prot = pd.DataFrame(df_prot.value_counts()).reset_index()
    df_prot.columns = ['Signature accession merged', 'Signature description merged', 'accession', 'organism_name', 'group', '#of proteins']
    df_prot = df_prot[['Signature accession merged', 'organism_name', 'Signature description merged', '#of proteins']]

    df_dom_prot = df_dom.merge(df_prot, how='left', on=['Signature accession merged', 'Signature description merged', 'organism_name'])

    return df_dom_prot
    


def reshape_interproscan_file(interpro_df):
    interpro_df['Start location'] = interpro_df['Start location'].astype(str) 
    interpro_df['Stop location'] = interpro_df['Stop location'].astype(str)

    interpro_df_pid = interpro_df.groupby('protein_id').agg({'Signature accession': lambda x: ';'.join(x)})
    interpro_df_pid = interpro_df_pid.reset_index()


    interpro_df_anno = interpro_df.groupby('protein_id').agg({'Signature description': lambda x: ';'.join(x)})
    interpro_df_anno = interpro_df_anno.reset_index()

    interpro_df = interpro_df_pid.join(interpro_df_anno.set_index('protein_id'), on='protein_id')
    #interpro_df = interpro_df[interpro_df['Signature accession'].str.contains(selected_domains[0])]
    ## one protein per line and the correct order, sort on protein and start values before joining the rows
    
    return interpro_df

def add_tm(df_matrix, tm_count):

    df_tm = pd.read_csv(str(tm_count), sep='\t')

    df_matrix = df_tm.merge(df_matrix, how='left', on='protein_id')
    df_matrix = df_matrix.fillna(0)

    return df_matrix

def merge_similar_domains(domains_to_merge, interpro_df):
    
    domains_to_merge_df = pd.read_csv(domains_to_merge, sep='\t')[['Signature accession','Signature accession merged', 'Signature description merged']]
    interpro_df_added = interpro_df.merge(domains_to_merge_df, how='left', on='Signature accession')
    interpro_df_added['Signature accession'] = interpro_df_added['Signature accession'].astype(str)

    return interpro_df_added


def main(interproscan_file, analysis, out_file_species, out_file_domains, merged_domains, protein_mapping, genome_stat, genome_stat_p_mapp_eu, domain_completness, tm_count, out_file_merged_domain_hmm_active_pkinase):
    
    # create an info file of genome info for pvc and eu
    df_info = create_df_info(protein_mapping, genome_stat, genome_stat_p_mapp_eu)
    
    # open interproscan and merge pfam domain ids based on input file
    interpro_original_df = open_interproscan_file(interproscan_file, analysis)
    interpro_original_df = interpro_original_df[['protein_id','Signature accession', 'Signature description']]
    interpro_original_df = merge_similar_domains(merged_domains, interpro_original_df)

    interpro_df =  interpro_original_df[['protein_id','Signature accession merged', 'Signature description merged']]
    
    # domain completness add full and fractional pkinase domain from the hmm search instead of the inerproscan pkinase domains
    #domain_completness_df = pd.read_csv(domain_completness, sep='\t')[['protein_id','domain_completness']]
    #domain_completness_df.columns = ['protein_id', 'Signature description merged']
    #domain_completness_df['Signature accession merged']='PF00069'    
    #interpro_df = interpro_df[interpro_df['Signature accession merged']!='PF00069']
    #interpro_df = pd.concat([interpro_df, domain_completness_df])
    
    # save active pkinases from hmm and their associated pfam domains with merged domain info
    #print(df_info.columns)
    interpro_df_info_added  = interpro_df.merge(df_info[['protein_id','accession','group', 'organism_name']], how='left',  on='protein_id')
    #i_full = interpro_df_info_added[interpro_df_info_added['Signature description merged']=='full'].index
    #i_fraciton = interpro_df_info_added[interpro_df_info_added['Signature description merged']=='full merged'].index
    #interpro_df_info_added.at[i_full,'Signature accession merged']='PF00069_full'
    #interpro_df_info_added.at[i_fraciton,'Signature accession merged']='PF00069_full'
    interpro_df_info_added.to_csv(out_file_merged_domain_hmm_active_pkinase, sep='\t', index= False)
    
    # count the domains and tm per protein
    df_t = pd.DataFrame(interpro_df.value_counts()).reset_index()
    df_t.columns = ['protein_id','Signature accession merged', 'Signature description merged','domain_count']  
    df = df_t[['protein_id','Signature description merged', 'domain_count']]
    df = df.pivot(index='protein_id',columns='Signature description merged', values='domain_count')
    df = df.fillna(0)
    df = df.merge(df_info[['protein_id','accession','group']], how='left', on='protein_id')
    df = add_tm(df, tm_count)
    f2 = df.to_csv(out_file_domains, sep='\t', index= False)
    
    # count the number of domains and number of proteins with these domains per species
    df_species_domain_count = create_dom_count(interpro_df, df_info)
    
    f1 = df_species_domain_count.to_csv(out_file_species, sep='\t', index= False)

    return f1, f2


# get input and run the script
parser = argparse.ArgumentParser()
parser.add_argument("--interpro_path")
parser.add_argument("--domain_completness")
parser.add_argument("--out_file_path_species")
parser.add_argument("--out_file_path_domains")
parser.add_argument("--analysis")
parser.add_argument("--merged_domains")
parser.add_argument("--protein_mapping")
parser.add_argument("--genome_stat")
parser.add_argument("--genome_stat_p_mapp_eu")
parser.add_argument("--tm_count")
parser.add_argument("--out_file_merged_domain_hmm_active_pkinase")

args = parser.parse_args()

#global selected_domains
#selected_domains = ['PF00069'] #['PF04542', 'PF08281', 'PF04539', 'PF04545']

main(args.interpro_path, args.analysis, args.out_file_path_species, args.out_file_path_domains, args.merged_domains, args.protein_mapping, args.genome_stat, args.genome_stat_p_mapp_eu, args.domain_completness, args.tm_count, args.out_file_merged_domain_hmm_active_pkinase)