import pandas as pd

df = pd.read_csv('../../results/02_orthofinder/single_copy_orthogroups.tsv', sep='\t', names=['n_species', 'orthogroup'])
orthogroups = df['orthogroup'].to_list()

rule all:
    input:
        expand("../../results/07_align_sco/trimmed_alignments/{orthogroup}.trimal.aln", orthogroup=orthogroups)

rule msa:
    """
    Align single-copy orthogroups
    """
    input:
        "../../results/02_orthofinder/orthofinder_results/Orthogroup_Sequences/{orthogroup}.fa"
    output:
        "../../results/07_align_sco/alignments/{orthogroup}.aln"
    threads:
        4
    conda:
        "../envs/mafft.yaml"
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule trim_alignment:
    """
    Filter the alignments using trimAl
    """
    input:
        "../../results/07_align_sco/alignments/{orthogroup}.aln"
    output:
        "../../results/07_align_sco/trimmed_alignments/{orthogroup}.trimal.aln"
    conda:
        "../envs/trimal.yaml"
    threads:
        1
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"
