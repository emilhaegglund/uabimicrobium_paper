import pandas as pd

# configuration
configfile: "config.yaml"

# add the input functions for this snakemake script
include: 'input-functions.smk'

# input list/file to run rules
pvc_df = pd.read_csv("../../data/accession_numbers.csv")
accession_list_pvc = pvc_df["Assembly Accession"].unique().tolist()
print(accession_list_pvc)
accession_list_eu = config["euk_accessions"]
print(accession_list_eu)
accession_all_list = accession_list_pvc + accession_list_eu
print("test")
groups = pvc_df['Group'].unique().tolist() + ['Nematoda', 'Ascomycota']

# wildcard constrainst
#accession_all_list_const = '|'.join(accession_all_list)
#group_const = '|'.join(groups)

def merge_hmm_folder_input(wildcards):
    input_files = []
    for accession in accession_list_eu:
        path = os.path.join("../../results/12-pkinase-analysis/hmmsearch-euk/", accession + "-pkinase-reformat.tsv")
        input_files.append(path)
    for accession in accession_list_pvc:
        path = os.path.join("../../results/12-pkinase-analysis/hmmsearch-pvc/", accession + "-pkinase-reformat.tsv")
        input_files.append(path)
    return input_files


rule all:
    input:
        '../../results/12-pkinase-analysis/active-sites/domain-completness-table.tsv',
        #'results/table/protein_kinase_table_merged.tsv',
        #'results/dag_figures/03_active_sites_filter.png'

rule merge_hmm_folder:
    input:
        merge_hmm_folder_input
    output:
        dir('../../results/12-pkinase-analysis/hmmsearch/')
    params:
        dir_path='../../results/12-pkinase-analysis/hmmsearch/'
    shell:
        """
        mkdir {params.dir_path}
        cp {input} {params.dir_path}
        """

rule domain_completness_table:
    input:
        hmm_folder = '../../results/12-pkinase-analysis/hmmsearch'
    output:
        output_1 = '../../results/12-pkinase-analysis/active-sites/domain-completness-table.tsv',
        output_2 = '../../results/12-pkinase-analysis/active-sites/domain-completness-table-active-and-not.tsv'
    params:
        group_list = groups,
        msa_folder = '../../results/12-pkinase-analysis/active-sites/domain-msa',
    script:
        'scripts/domain_type.py'

#rule pkinase_merged_table:
#    input:
#        'results/active_sites/domain_completness_table.tsv'
#    output:
#        'results/active_sites/domain_completness_table_merged.tsv'
#    script:
#        'scripts/merged_pkinase_domains.py'
#
#rule create_eu_genomes_table:
#    params:
#        'results/ref_proteomes/'
#    output:
#        'results/table/protein_mapping_other_species.tsv'
#    script:
#        'scripts/create_eu_genomes_table.py'
#
#rule pkinase_abundance_table_with_full:
#    input:
#        p_count = config['p_count'],
#        genome_stat = config['genome_statistics'],
#        genome_stat_eu = 'results/table/protein_mapping_other_species.tsv',
#        OG_orthomcl = config['ogs'],
#        OG_orthofinder = config['ogs_orthofinder'],
#        OG_orthofinder_unassigned = config['unassigned_orthofinder'],
#        p_mapp = config['protein_mapping'],
#        pkinase_active = 'results/active_sites/domain_completness_table_merged.tsv',
#        fasta_file = config['merged_pvc_and_eu_pkinases']
#    output:
#        output = 'results/table/protein_kinase_table_merged.tsv'
#    params:
#        hmm_path = 'results/hmm/'
#    script:
#        'scripts/pkinase_species_table.py'
#
#rule create_dag_figure:
#    output:
#        'results/dag_figures/03_active_sites_filter.png'
#    shell:
#        'snakemake -s snakefiles/03_active_sites_filter.smk --forceall --rulegraph | dot -Tpng > {output}'
