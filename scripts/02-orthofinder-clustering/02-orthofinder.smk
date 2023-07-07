# Run orthofinder on the proteomes in the ../results/01_download_data folder
# Also extract the single protein copy orthogroups from the orthofinder results, at least 95% of the genomes must have a protein in the OG
rule all:
    input:
        '../../results/02_orthofinder/orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv',
        '../../results/02_orthofinder/single_copy_orthogroups.tsv'

rule orthofinder:
    """
    Run OrthoFinder
    """
    input:
        '../../results/01_download_data/genbank_proteomes/'
    output:
        '../../results/02_orthofinder/orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv'
    params:
        outdir='../../results/02_orthofinder/orthofinder_results'
    conda:
        '../envs/orthofinder.yaml'
    threads:
        28
    shell:
        """
        rm -rf {params.outdir};
        rm -rf {input}/OrthoFinder
        orthofinder -os -M msa -S diamond -t 28 -a 28 -f {input};
        mv {input}/OrthoFinder/Results_* {params.outdir};
        rm -rf {input}/OrthoFinder
        """

rule get_single_copy:
    input:
        '../../results/02_orthofinder/orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv'
    output:
        '../../results/02_orthofinder/single_copy_orthogroups.tsv'
    threads:
        1
    shell:
        'python orthogroup-selection.py {input} {output}'
