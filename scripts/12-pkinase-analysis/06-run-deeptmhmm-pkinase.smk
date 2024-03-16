"""
Snakemake script to run phobius locally on the pkinases of pvc and other species (e.coli and yeast/c elegans.
"""


import pandas as pd

# configuration
configfile: '02_config.yml' 
    
rule all:
    input:
        'results/phobius/tm_rd_barplot_table.tsv',
        'results/dag_figures/06_run_phobius_pkinase.png'
        
rule run_phobius:
    input:
        'results/interproscan/merged_pvc_and_eu_pkinases.fa'
    output:
        'results/phobius/merged_pvc_and_eu_pkinases.fa.tsv'
    log:
        'logs/phobius.log'
    shell:
        'phobius.pl {input} > {output} 2> {log}'
        
rule generate_phobius_tables:
    input:
        phobius_file = 'results/phobius/merged_pvc_and_eu_pkinases.fa.tsv',
        genome_stat = config['genome_statistics'],
        protein_kinase_table = 'results/table/protein_kinase_table_merged.tsv'
    output:
        phobius_reformat = 'results/phobius/phobius_reformatted.tsv',
        tm_protein_table = 'results/phobius/TM_comp.tsv',
        tm_species_table = 'results/phobius/pkinase_phobius.tsv'
    script:
        'scripts/tm_abundance_table.py'

rule generate_kdd_rd_table:
    input:
        phobius = 'results/phobius/pkinase_phobius.tsv',
        active_sites_table = 'results/active_sites/domain_completness_table_merged.tsv',
        p_mapp = config['protein_mapping'],
        genome_stat = config['genome_statistics'],
        genome_stat_eu = config['eu_genomes']
    output:
        kdd_rd_table = 'results/phobius/tm_rd_barplot_table.tsv',
        tm_figure_name = 'results/phobius/tm_rd_barplot.png'
    script:
        'scripts/KDD_RD.py'

rule create_dag_figure:
    output:
        'results/dag_figures/06_run_phobius_pkinase.png'
    shell:
        'snakemake -s snakefiles/06_run_phobius_pkinase.smk -F --rulegraph | dot -Tpng > {output}'
