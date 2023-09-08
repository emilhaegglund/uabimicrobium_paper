import pandas as pd

pvc_df = pd.read_csv("../../data/accession_numbers.csv", sep=",")
pvc_accessions = pvc_df["Assembly Accession"].unique()
pvc_groups = pvc_df["Group"].unique()

euk_accessions = ["GCA_000002985.3",
                  "GCA_000002945.2",
                  "GCA_000146045.2"]

def merge_fasta_input(wildcards):
    """ Function to create input files for each group in the merge step """
    group_accessions = pvc_df[pvc_df["Group"] == wildcards.group]["Assembly Accession"].unique()
    input_files = []
    for accession in group_accessions:
        infile = os.path.join("../../results/12-pkinase-analysis/active-sites/domain-fasta-pvc/", accession + "-PF00069.faa")
        input_files.append(infile)
    return input_files

rule all:
    input:
        expand('../../results/12-pkinase-analysis/active-sites/domain-msa/PF00069-{group}.aln', group=pvc_groups),
        '../../results/12-pkinase-analysis/active-sites/domain-msa/PF00069-euk.aln',
        "../../snakemake-dags/02-activet-sites.png"


rule create_domain_fasta_files_pvc:
    input:
        hmm_result = "../../results/12-pkinase-analysis/hmmsearch-pvc/{accession}-pkinase.tsv",
        fasta = "../../results/01_download_data/genbank_proteomes/{accession}.faa"
    output:
        '../../results/12-pkinase-analysis/active-sites/domain-fasta-pvc/{accession}-PF00069.faa'
    shell:
        'python scripts/hmm_create_domain_seq_files.py --hmm_domain_info_file {input.hmm_result} --species_fasta_file {input.fasta} --species_domain_out_fasta_file {output}'


rule merge_fasta_on_group_pvc:
    input:
        merge_fasta_input
    output:
        '../../results/12-pkinase-analysis/active-sites/domain-fasta-groups/PF00069-{group}.faa'
    shell:
        'cat {input} > {output}'

rule create_domain_fasta_files_euk:
    input:
        hmm_result = "../../results/12-pkinase-analysis/hmmsearch-euk/{accession}-pkinase.tsv",
        fasta = "../../results/12-pkinase-analysis/euk-proteomes/{accession}.faa"
    output:
        '../../results/12-pkinase-analysis/active-sites/domain-fasta-euk/{accession}-PF00069.faa'
    shell:
        'python scripts/hmm_create_domain_seq_files.py --hmm_domain_info_file {input.hmm_result} --species_fasta_file {input.fasta} --species_domain_out_fasta_file {output}'

rule merge_fasta_on_group_euk:
    input:
        expand("../../results/12-pkinase-analysis/active-sites/domain-fasta-euk/{accession}-PF00069.faa", accession=euk_accessions)
    output:
        '../../results/12-pkinase-analysis/active-sites/domain-fasta-groups/PF00069-euk.faa'
    shell:
        'cat {input} > {output}'

rule msa_with_mafft:
    input:
        '../../results/12-pkinase-analysis/active-sites/domain-fasta-groups/PF00069-{group}.faa'
    output:
        '../../results/12-pkinase-analysis/active-sites/domain-msa/PF00069-{group}.aln'
    conda:
        '../envs/mafft.yaml'
    log:
        '../../logs/02-active-sites/mafft/domain-merging-{group}.log'
    threads:
        12
    shell:
        'mafft-linsi --thread {threads} {input} > {output} 2> {log}'

rule create_dag_figure:
    output:
        "../../snakemake-dags/02-activet-sites.png"
    shell:
        'snakemake -s 02-active-sites.smk --forceall --rulegraph | dot -Tpng > {output}'
