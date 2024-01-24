"""
Author: Emil Hagglund
Date: 240124

Workflow to perform phylogenetic analysis of the Pkinase domain inside the
PVC.
"""
rule all:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.bmge_blosum30.aln"

rule extract_pkinase_seq:
    """Use inforrmation from Supplementary Table 5 to extract the protein
    sequences for the Pkinase domains from PVC genomes and write them to a
    fasta file."""
    input:
        suppl_table_5 = "../../tables/supplementary-table-5.tsv",
        proteoeme_dir = "../../results/01_download_data/genbank_proteomes/"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.faa"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python extract-pkinase-sequences.py {input.suppl_table_5} {input.proteoeme_dir} {output}"

rule align_pkinase:
    """Use Mafft-linsi to align all Pkinase sequences."""
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.faa"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.aln"
    conda:
        "../envs/mafft.yaml"
    threads:
        12
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule bmge_trimming:
    """Trim the alignment using BMGE."""
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.bmge_blosum30.aln"
    conda:
        "../envs/bmge.yaml"
    shell:
        "bmge -i {input} -of {output} -t AA -m BLOSUM30"

rule fasttree:
    """First create an initial phylogeny using FastTree."""
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.bmge_blosum30.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.bmge_blosum30.fasttree.nwk"
    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"