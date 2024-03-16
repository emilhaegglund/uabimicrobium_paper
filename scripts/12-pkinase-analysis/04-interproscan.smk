
import pandas as pd


# configuration
configfile: '02_config.yml'

accession_list_pvc = pd.read_csv(config['p_count'], sep='\t')['assembly_accession'].values.tolist()
accession_list_eu = config['eu_accession']
accession_list_all = accession_list_pvc + accession_list_eu

rule all:
    input:
        'results/dag_figures/04_interproscan.png',
        'results/interproscan/merged_pvc_and_eu_pkinases_filtered.tsv'

rule extract_pvc_pkinase_aaseq:
    input:
        'results/hmm/{accession}_pkinase_reformated.tsv',
        config['proteomes_folder'] + '{accession}.faa'
    output:
        'results/interproscan/{accession}_pkinases.fa'
    wildcard_constraints:
        accession = '|'.join(accession_list_pvc)
    shell:
        'python scripts/interpro_merge_fasta.py --hmm_file {input[0]} --out_file {output[0]} --proteome_path {input[1]}'

rule extract_eu_pkinase_aaseq:
    input:
        'results/hmm/{accession}_pkinase_reformated.tsv',
        config['ref_proteom_folder'] + '{accession}.faa'
    output:
        'results/interproscan/{accession}_pkinases.fa'
    wildcard_constraints:
        accession = 'GCA_000005845.2|GCA_000002985.3|GCA_000146045.2|GCA_000002945.2'
    shell:
        'python scripts/interpro_merge_fasta.py --hmm_file {input[0]} --out_file {output[0]} --proteome_path {input[1]}'

rule merge_eu_pvc_fasta:
    input:
        expand('results/interproscan/{accession}_pkinases.fa', accession = accession_list_all)
    output:
        'results/interproscan/merged_pvc_and_eu_pkinases.fa'
    shell:
        'cat {input} >> {output}'

rule run_interproscan:
    input:
        'results/interproscan/merged_pvc_and_eu_pkinases.fa'
    output:
        'results/interproscan/merged_pvc_and_eu_pkinases.tsv'
    conda:
        '../env/interproscan.yml'
    log:
        'logs/interproscan_pkinase.log'
    shell:
        'interproscan.sh -appl CDD,COILS,Pfam,SMART,SUPERFAMILY,TIGRFAM -i {input} -f TSV -o {output} >> {log}'

rule filter_interproscan:
    input:
        interproscan_raw = 'results/interproscan/merged_pvc_and_eu_pkinases.tsv',
        pkinase_table = 'results/active_sites/domain_completness_table_merged.tsv'
    output:
        filtered = 'results/interproscan/merged_pvc_and_eu_pkinases_filtered.tsv'
    script:
        'scripts/interproscan_filtered.py'

rule create_dag_figure:
    output:
        'results/dag_figures/04_interproscan.png'
    shell:
        'snakemake -s snakefiles/04_interproscan.smk --forceall --rulegraph | dot -Tpng > {output}'

