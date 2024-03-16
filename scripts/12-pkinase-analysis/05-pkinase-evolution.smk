import pandas as pd

# configuration
configfile: '02_config.yml'

# add the input functions for this snakemake script
include: 'input_functions.smk'

domain_info_df = pd.read_csv('results/active_sites/domain_completness_table_merged.tsv', sep='\t')
groups = domain_info_df['group'].drop_duplicates().values.tolist()

rule all:
    input:
        'results/gain_loss_protein_kinase_tree.png',
        'results/dag_figures/05_pkinase_evolution.png'

rule protein_kinase_gain_loss_tree:
    input:
        tree = config['count_tree'],
        pkinase_species_result = 'results/table/protein_kinase_table_merged.tsv',
        gained_og = config['gained_og'],
        lost_og = config['lost_og'],
        genome_stat = config['genome_statistics']
    output:
        out_file_name = 'results/gain_loss_protein_kinase_tree.png'
    script:
        'scripts/extract_gains_loss_posterior.py'

rule create_dag_figure:
    output:
        'results/dag_figures/05_pkinase_evolution.png'
    shell:
        'snakemake -s snakefiles/05_pkinase_evolution.smk --forceall --rulegraph | dot -Tpng > {output}'
