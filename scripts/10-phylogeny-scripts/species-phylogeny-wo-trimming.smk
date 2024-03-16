
replicates = list(range(100))
rule all:
    input:
        # Concatenated alignment
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.aln",
        # Guide tree for PMSF
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.treefile",
        # ModelFinder for MixtureModel
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.modelfinder.log",
        # ML-tree using PMSF of MixtureModel
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.ml.treefile",
        # Bootstrap trees (run this independently to speed up calculations)
        expand("../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.boot{replicate}.boottrees", replicate=replicates),
        # ML-tree with mapped support values
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.ml-boot.suptree"

rule concatenate_sco:
    input:
        protein_mappings="../../results/01_download_data/protein_mapping.tsv",
        alignment_dir="../../results/06_align_sco/alignments/",
        orthogroups="../../results/07_species_phylogeny_single_copy_orthogroups/orthogroups.txt"
    output:
        concat="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.aln"
    shell:
        """
        python ../concat_alignment.py --protein_mapping {input.protein_mappings} \
            --orthogroups {input.orthogroups} \
            --alignments {input.alignment_dir} \
            --concatenated {output.concat}
        """

rule species_phylogeny_guide_tree:
    """
    As a first step create a guide tree to use for the PMSF approximation.
    """
    input:
        concat="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.aln",
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.treefile"
    params:
        pre="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming",
        uf=1000,
        alrt=1000
    conda:
        "../envs/iqtree.yaml"
    threads:
        30
    shell:
        "iqtree -s {input.concat} -m LG+F+G -seed 123 -pre {params.pre} -bnni -bb {params.uf} -alrt {params.alrt} -nt {threads} -redo"

rule species_phylogeny_model_finder:
    """
    Run ModelFinder on the LG+C10-C60 series to find out what model to used in the downstream steps.
    Results:
    """
    input:
        concat="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.aln",
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.modelfinder.log"
    params:
        pre="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.modelfinder"
    conda:
        "../envs/iqtree.yaml"
    threads:
        30
    shell:
        """
        iqtree -s {input.concat} \
            -m TESTONLY \
            -mset LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 \
            -pre {params.pre} \
            -seed 123 \
            -nt {threads} \
            -redo
        """

rule species_phylogeny_pmsf_site_freq:
    """
    Infer site-specific frequency profiles.
    """
    input:
        concat="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.aln",
        guide_tree="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.treefile"
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.sf.sitefreq"
    params:
        pre="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.sf",
        model="LG+C60+F+G"
    conda:
        "../envs/iqtree.yaml"
    threads:
        30
    shell:
        """
        iqtree -s {input.concat} \
            -m {params.model} \
            -pre {params.pre} \
            -ft {input.guide_tree} \
            -n 0 \
            -seed 123 \
            -nt {threads} \
            -redo
        """

rule species_phylogeny_pmsf_ml:
    """
    Infer the ML-phylogeny using the site-frequencies and substitution model selected in the previous steps.
    """
    input:
        concat="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.aln",
        site_frequency="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.sf.sitefreq"
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.ml.treefile"
    params:
        pre="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.ml",
        model="LG+C60+F+G"
    conda:
        "../envs/iqtree.yaml"
    threads:
        6
    shell:
        """
        iqtree -s {input.concat} \
            -m {params.model} \
            -pre {params.pre} \
            -fs {input.site_frequency} \
            -nt {threads} \
            -seed 123 \
            -redo
        """

rule species_phylogeny_pmsf_bootstrap:
    """
    Try to speed up the computation time by separating the bootstrap tree on different computers.
    Skip setting seed, don't know if it will generate same bootstrap alignment with a set seed.
    """
    input:
        concat="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.aln",
        site_frequency="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.sf.sitefreq"
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.boot{replicate}.boottrees"
    params:
        pre="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.boot{replicate}",
        model="LG+C60+F+G"
    conda:
        "../envs/iqtree.yaml"
    threads:
        6
    shell:
        """
        iqtree -s {input.concat} \
            -m {params.model} \
            -pre {params.pre} \
            -fs {input.site_frequency} \
            -bo 1 \
            -nt {threads} \
            -redo
        """

rule merge_boottrees:
    input:
        expand("../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.boot{replicate}.boottrees", replicate=replicates)
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/alltrees"
    shell:
        "cat {input} > {output}"

rule map_support_values:
    input:
        ml_tree = "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.ml.treefile",
        boottrees = "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/alltrees"
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.ml-boot.suptree"
    params:
        pre="../../results/07_species_phylogeny_single_copy_orthogroups/no_trimming/species_phylogeny.no_trimming.ml-boot",
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        iqtree -sup {input.ml_tree} -t {input.boottrees} -pre {params.pre}
        """
