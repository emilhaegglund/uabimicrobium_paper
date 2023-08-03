#dataset = ['11_species_phylogeny_aminoGC/aminoGC_high']
dataset = ['07_species_phylogeny_single_copy_orthogroups']
replicates = list(range(100))

rule all:
    input:
        expand("../../results/{dataset}/bmge/species_phylogeny.treefile", dataset=dataset),
        expand("../../results/{dataset}/bmge/species_phylogeny.C60.sitefreq", dataset=dataset),
        expand("../../results/{dataset}/bmge/species_phylogeny.C60.boot{replicate}.boottrees", dataset=dataset, replicate=replicates),
        expand("../../results/{dataset}/bmge/species_phylogeny.C60.alltrees_blength.contree.treefile", dataset=dataset)

rule concatenate_sco:
    input:
        protein_mappings="../../results/01_download_data/protein_mapping.tsv",
        alignment_dir="../../results/07_align_sco/alignments",
        orthogroups="../../results/{dataset}/orthogroups.txt"
    output:
        concat="../../results/{dataset}/bmge/species_phylogeny.aln"
    shell:
        "python ../concat_alignment.py --protein_mapping {input.protein_mappings} --orthogroups {input.orthogroups} --alignments {input.alignment_dir} --concatenated {output.concat}"

rule bmge_stationary_trimming:
    input:
        "../../results/{dataset}/bmge/species_phylogeny.aln"
    output:
        "../../results/{dataset}/bmge/species_phylogeny.trimmed.aln"
    conda:
        "../envs/bmge.yaml"
    shell:
        "bmge -i {input} -of {output} -t AA -s FAST -h 0:1 -g 1"

rule species_phylogeny:
    input:
        concat="../../results/{dataset}/bmge/species_phylogeny.trimmed.aln",
    output:
        "../../results/{dataset}/bmge/species_phylogeny.treefile"
    params:
        pre="../../results/{dataset}/bmge/species_phylogeny",
        uf=1000,
        alrt=1000
    conda:
        "../envs/iqtree.yaml"
    threads:
        8
    shell:
        "iqtree -s {input.concat} -m LG+F+G -pre {params.pre} -bnni -bb {params.uf} -alrt {params.alrt} -nt {threads} -redo"

# very heavy, need 79517 MB RAM. You might want to do this step on a high-memory server
rule species_phylogeny_pmsf_pre_step:
    input:
        concat="../../results/{dataset}/bmge/species_phylogeny.trimmed.aln",
        guide_tree="../../results/{dataset}/bmge/species_phylogeny.treefile"
    output:
        "../../results/{dataset}/bmge/species_phylogeny.C60.sitefreq"
    params:
        pre="../../results/{dataset}/bmge/species_phylogeny.C60"
    conda:
        "../envs/iqtree.yaml"
    threads:
        8
    shell:
        "iqtree -s {input.concat} -m LG+C60+F+G -ft {input.guide_tree} -pre {params.pre} -n 0 -nt {threads} -redo"

# Use the sitefreq file from the step above
rule species_phylogeny_pmsf:
    input:
        concat="../../results/{dataset}/bmge/species_phylogeny.trimmed.aln",
        site_frequency="../../results/{dataset}/bmge/species_phylogeny.C60.sitefreq"
    output:
        "../../results/{dataset}/bmge/species_phylogeny.C60.boot{replicate}.boottrees"
    params:
        pre="../../results/{dataset}/bmge/species_phylogeny.C60.boot{replicate}"
    conda:
        "../envs/iqtree.yaml"
    threads:
        6
    shell:
        "iqtree -s {input.concat} -m LG+C60+F+G -pre {params.pre} -fs {input.site_frequency} -bo 1 -nt {threads} -redo"

rule consensus_tree_step1:
    input:
        expand("../../results/{dataset}/bmge/species_phylogeny.C60.boot{replicate}.boottrees", replicate=replicates, dataset=dataset)
    output:
        "../../results/{dataset}/bmge/species_phylogeny.C60.alltrees"
    shell:
        "cat ../../results/{dataset}/bmge/species_phylogeny.C60.boot*.boottrees > {output}"

rule consensus_tree_step2:
    input:
        "../../results/{dataset}/bmge/species_phylogeny.C60.alltrees"
    output:
        "../../results/{dataset}/bmge/species_phylogeny.C60.alltrees.contree"
    conda:
        "../envs/iqtree.yaml"
    shell:
        "iqtree -con -t {input}"

rule consensus_branch_length:
    input:
        alignment = "../../results/{dataset}/bmge/species_phylogeny.trimmed.aln",
        consensus_tree = "../../results/{dataset}/bmge/species_phylogeny.C60.alltrees.contree",
    params:
        pre = "../../results/{dataset}/bmge/species_phylogeny.C60.alltrees_blength.contree"
    output:
        "../../results/{dataset}/bmge/species_phylogeny.C60.alltrees_blength.contree.treefile"
    conda:
        "../envs/iqtree.yaml"
    threads:
        8
    shell:
        "iqtree -s {input.alignment} -te {input.consensus_tree} -pre {params.pre} -m LG -nt {threads}"
