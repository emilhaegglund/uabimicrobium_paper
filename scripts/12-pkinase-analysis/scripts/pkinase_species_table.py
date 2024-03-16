import pandas as pd
import os
from Bio import SeqIO


def read_hmm_file(hmm_file, hmm_path):
    col_names = ['target name','accession_1', 'query name', 'accession_2', 'E-value_1', 'score_1', 'bias_1', 'E-value_2', 'score_2', 'bias_2', 'exp', 'reg', 'clu', 'ov', 'env', 'dom','rep', 'inc' ,'description of target']
    df = pd.read_csv(hmm_path + hmm_file, names = col_names, sep='\t', comment='#')
    return df

def open_file(my_file):
    """
    open the hmm file as a dataframe
    """

    col_names = ['target name','accession','tlen', 'query name','accession domain','qlen','E-value', 'score', 'bias','#', 'of', 'c-Evalue', 'i-Evalue', 'score2', 'bias2', 'domain from', 'domain to', 'protein from', 'protein to']
    
    df= pd.read_csv(my_file, sep='\s',comment='#', header=None, names= col_names, usecols=col_names, index_col=False)  
    df = df.fillna('Missing')

    return df

def read_interpro(interpro_path):
    col_names = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession', 'Signature description', 'Start location', 'Stop location', 'Score', 'Status',
        'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
        'optional: GO terms', 'optional: Pathway annotation']
    all_interproscan_raw_df = pd.read_csv(interpro_path, sep='\t', header=None, names=col_names)
    return all_interproscan_raw_df


def read_og_file(og_file_path, clustering_software):
    og_df = pd.read_csv(og_file_path, sep='\s', header=None)
    og_df = og_df.melt(id_vars=0).reset_index(drop=True)
    og_df = og_df.dropna() # remove lines with missing values after the melt command
    og = og_df[[0, 'value']] # extract the orthogroup and protein columns
    og.columns = [ str(clustering_software) + '_orthogroup', 'sseqid'] # rename the orthogroup and protein columns (sseqid will be used as key for mergin to blast results)
    return og

def read_og_unassigned_file(og_file_path, clustering_software):

    og_df = pd.read_csv(og_file_path, sep='\t', keep_default_na=False)
    og = og_df['Orthogroup'].values.tolist()
    #print(og)
    #og_df = og_df.melt(id_vars=0).reset_index(drop=True)
    #og_df = og_df.dropna() # remove lines with missing values after the melt command
    #og = og_df[[0, 'value']] # extract the orthogroup and protein columns
    #og.columns = [ str(clustering_software) + '_orthogroup', 'sseqid'] # rename the orthogroup and protein columns (sseqid will be used as key for mergin to blast results)
    return og

def get_ogs(og_df, protein_ids, clustering_software, orthofinder_unassigned):

    og_list = og_df[og_df['sseqid'].isin(protein_ids)][str(clustering_software) + '_orthogroup'].drop_duplicates().values.tolist()
    og_proteins = og_df['sseqid'].values.tolist()
    #og_proteins_orthofinder_un = orthofinder_unassigned_df['sseqid'].values.tolist()
    singletons = []
    
    if clustering_software == 'orthomcl':
        for protein in protein_ids:
            if  protein.rstrip() not in og_proteins:
                singletons.append(protein)

    if clustering_software == 'orthofinder':
        for OG in og_list:
            if OG.rstrip().split(':')[0] in orthofinder_unassigned:
                protein = og_df[og_df[str(clustering_software) + '_orthogroup']==OG]['sseqid']
                singletons.append(protein)

    return og_list, singletons

def get_annotation(fasta_file, protein_ids):
    protein_ids = protein_ids.split(';')
    annotation_temp = []
    current_anno = ''

    fasta_records = open_fasta(fasta_file)
    for record in fasta_records:
        if record.id in protein_ids: 
                current_anno = record.description.split(record.id)[-1].split('[')[0].rstrip().lstrip()
                annotation_temp.append([record.id, current_anno])
    
    annotation_df = pd.DataFrame(annotation_temp)
    annotation_df.columns = ['protein_id', 'annotation']
    
    annotation = []
    for pid in protein_ids:
        annotation.append(annotation_df[annotation_df['protein_id']==pid]['annotation'].values.tolist()[0])
    annotation = ';'.join(annotation)
    return annotation

def open_fasta(fasta_file):
    records = SeqIO.parse(fasta_file, 'fasta')
    return records

def get_protein_count(fasta_file):

    number_of_proteins = 0
    fasta_records = open_fasta(fasta_file)

    for record in fasta_records:
        number_of_proteins +=1

    return number_of_proteins

def get_id_species_table(df_pvc, protein_mapping):
    pkinase_pids = df_pvc[df_pvc['Analysis']=='Pfam']
    pkinase_pids = pkinase_pids[pkinase_pids['Signature accession']=='PF00069']['protein_id'].drop_duplicates().values.tolist()
    
    pkinase_domain_table = df_pvc[df_pvc['protein_id'].isin(pkinase_pids)]
    #add species to the pkinase_domain_table
    pid_species_df = pkinase_domain_table.merge(protein_mapping, how='left', on='protein_id')
    return pid_species_df

def generate_pkinase_table(genome_acc, genome_stat_df,hmm_genome_stat_df, og_df_orthomcl, og_df_orthofinder, orthofinder_unassigned_df, fasta_file):
    
    pkinase_table = []
    

    # create a table with species and their protein ids
    species_acc_list = hmm_genome_stat_df['accession'].drop_duplicates().values.tolist()

    #for species_acc in species_acc_list:
    for species_acc in genome_acc:
        
        protein_count = genome_stat_df[genome_stat_df['accession']==species_acc]['protein count'].values[0]
        species = genome_stat_df[genome_stat_df['accession']==species_acc]['organism name'].values[0]

        if species_acc == 'GCA_009002475.1':
            group = 'Orphan_Uab'
        else:
            group = genome_stat_df[genome_stat_df['accession']==species_acc]['group'].values[0]
                
        if species_acc in species_acc_list:
            protein_ids = hmm_genome_stat_df[hmm_genome_stat_df['accession']==species_acc]['protein_id'].drop_duplicates().values.tolist()
            protein_ids_reformat = hmm_genome_stat_df[hmm_genome_stat_df['accession']==species_acc]['protein_id'].drop_duplicates().sort_values().str.cat(sep=';').rstrip()

            OGs_orthomcl, singletons_orthomcl = get_ogs(og_df_orthomcl, protein_ids, 'orthomcl', orthofinder_unassigned_df)
            OG_orthomcl_display = og_df_orthomcl[og_df_orthomcl['sseqid'].isin(protein_ids)]['orthomcl_orthogroup'].drop_duplicates().sort_values().str.cat(sep=';').rstrip()

            OGs_orthofinder, singletons_orthofinder = get_ogs(og_df_orthofinder, protein_ids, 'orthofinder', orthofinder_unassigned_df)

            annotation = get_annotation(fasta_file, protein_ids_reformat)
            frac = round((len(protein_ids)/protein_count)*100)
            # Add info to the pkinase table:
            pkinase_table.append([species_acc, species, group, len(protein_ids), protein_count, frac, len(OGs_orthomcl), len(singletons_orthomcl), len(OGs_orthofinder), len(singletons_orthofinder), annotation, protein_ids_reformat, OG_orthomcl_display])
            
        else:
            pkinase_table.append([species_acc, species, group, 0, protein_count, 0, 0, 0, 0, 0, '-', '-', '-'])  

    return pkinase_table

def save_table(pkinase_table, out_path):
    pkinase_df = pd.DataFrame(pkinase_table, columns=['accession', 'species', 'group', 'nr of pkinases', 'all proteins', 'fraction of pkinases (%)', 'orthomcl OGs', 'number of orthomcl singletons', 'orthofinder OGs', 'number of orthofinder singletons', 'annotation', 'pkinase protein ids', 'orthomcl OGs'])
    f = pkinase_df.to_csv(out_path, index=False, sep='\t')
    return f

def merge_files(protein_count_file, genome_stat_pvc, p_mapp_pvc, genome_stat_eu):
    protein_count_df = pd.read_csv(protein_count_file, sep='\t')[['assembly_accession', 'nr_proteins']]
    genome_stat_df_pvc = pd.read_csv(genome_stat_pvc, sep='\t')[['assembly_accession', 'organism_name', 'group']]
    genome_stat_df_pvc = genome_stat_df_pvc.merge(protein_count_df, how='left', on='assembly_accession')
    genome_stat_df_pvc.columns = ['accession', 'organism name', 'group', 'protein count']
    p_mapp = pd.read_csv(p_mapp_pvc, sep='\t')[['protein_id', 'accession']]
    genome_stat_df_pvc = genome_stat_df_pvc.merge(p_mapp, how='left', on='accession')
    genome_stat_df_pvc = genome_stat_df_pvc[['accession', 'protein_id', 'group', 'organism name','protein count']]
    genome_stat_df_eu = pd.read_csv(genome_stat_eu, sep='\t')[['accession', 'protein_id', 'group', 'organism name','protein count']]
    
    genome_stat = pd.concat([genome_stat_df_pvc, genome_stat_df_eu], ignore_index=True)
    genome_acc = genome_stat['accession'].drop_duplicates().values.tolist()
    
    return genome_stat, genome_acc

def merge_hmm(hmm_path):
    hmm_count = 0
    for hmm_file in os.listdir(hmm_path):
        if 'reformated' in hmm_file:
            if hmm_count == 0:
                df_hmm = open_file(hmm_path + hmm_file)
                hmm_count+=1
            else:
                df = open_file(hmm_path + hmm_file)
                df_hmm = pd.concat([df_hmm, df])
    df_hmm = df_hmm[['target name', 'tlen', 'query name', 'accession domain',
       'qlen', '#', 'of', 'domain from', 'domain to', 'protein from',
       'protein to']]
    df_hmm.columns = ['protein_id', 'tlen', 'domain name', 'domain accession',
       'qlen', '#', 'of', 'domain from', 'domain to', 'protein from',
       'protein to']
    return df_hmm

def main(hmm_path, protein_count_file, genome_stat_pvc, genome_stat_eu, OG_orthomcl, OG_orthofinder, OG_orthofinder_unassigned, p_mapp_pvc, pkinase_active, fasta_file, output):

    # open input files and generate dataframes
    genome_stat_df, genome_acc = merge_files(protein_count_file, genome_stat_pvc, p_mapp_pvc, genome_stat_eu)

    hmm_df = merge_hmm(hmm_path)
    pid_active = pd.read_csv(pkinase_active, sep='\t')['protein_id'].values.tolist()
    hmm_df = hmm_df[hmm_df['protein_id'].isin(pid_active)]

    hmm_genome_stat_df = hmm_df.merge(genome_stat_df, how='left', on='protein_id')

    og_df_orthomcl = read_og_file(OG_orthomcl, 'orthomcl')
    og_df_orthofinder = read_og_file(OG_orthofinder, 'orthofinder')
    orthofinder_unassigned_df = read_og_unassigned_file(OG_orthofinder_unassigned, 'orthofinder')

    # generate the pkinase table/output file, extract the protein hits with the hmm PF000069 profile found by hmm
    # also add information about the hits and make one row per species
    print('here')
    pkinase_table = generate_pkinase_table(genome_acc, genome_stat_df, hmm_genome_stat_df, og_df_orthomcl, og_df_orthofinder, orthofinder_unassigned_df, fasta_file)
    
    t = save_table(pkinase_table, output)

    return t

### run the script:

#input
hmm_path = str(snakemake.params.hmm_path)
protein_count_file = str(snakemake.input.p_count)
genome_stat_pvc = str(snakemake.input.genome_stat) 
genome_stat_eu = str(snakemake.input.genome_stat_eu) # ../results/protein_mapping_other_species.tsv'?
OG_orthomcl = str(snakemake.input.OG_orthomcl)
OG_orthofinder = str(snakemake.input.OG_orthofinder)
OG_orthofinder_unassigned = str(snakemake.input.OG_orthofinder_unassigned)
p_mapp_pvc = str(snakemake.input.p_mapp)
pkinase_active = str(snakemake.input.pkinase_active)
fasta_file = str(snakemake.input.fasta_file)
output = str(snakemake.output.output)

# run the main function
table = main(hmm_path, protein_count_file, genome_stat_pvc, genome_stat_eu, OG_orthomcl, OG_orthofinder, OG_orthofinder_unassigned, p_mapp_pvc, pkinase_active, fasta_file, output)



