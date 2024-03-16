# input functions for the snakemake script 02
import pandas as pd
def get_proteome_ncbi_path(wildcards):
    file_path = ''
    col_names = ['assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid','species_taxid', 'organism_name', 'infraspecific_name', 'isolate', 'version_status', 'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter', 'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path', 'excluded_from_refseq', 'relation_to_type_material', 'asm_not_live_date']
    asm_file = pd.read_csv('results/assembly_summary_genbank.txt', sep='\t', header=None, names=col_names)
    file_path = str(asm_file[asm_file['assembly_accession']==str(wildcards)]['ftp_path'].values[0])
    file_path += '/' + file_path.split('/')[-1]
    return file_path

def get_fasta_path(wildcards):
    if str(wildcards) in accession_list_eu:
        return config['ref_proteom_folder'] + '{accession}.faa'
    if str(wildcards) in accession_list_pvc:
        return config['proteomes_folder'] + '{accession}.faa'
    
def merge_fasta_get_group(wildcards):
    group_files = []
    
    if str(wildcards).strip() in pd.read_csv(config['eu_genomes'], sep='\t')['group'].drop_duplicates().values.tolist():
        df_eu = pd.read_csv(config['eu_genomes'], sep='\t')
        accessions = df_eu[df_eu['group']==str(wildcards)]['accession'].drop_duplicates().values.tolist()
    else:
        df_pvc = pd.read_csv(config['genome_statistics'], sep='\t')
        accessions = df_pvc[df_pvc['group']==str(wildcards)]['assembly_accession'].drop_duplicates().values.tolist()
    
    for acc in accessions:     
        group_files.append('results/active_sites/domain_fasta/PF00069_' + str(acc) + '.faa')
    
    return group_files

