orthogroups = glob_wildcards('../../results/06_align_sco/alignments/{orthogroup}.aln').orthogroup

rule all:
    input:
        "../../results/09_species_phylogeny_aminoGC/aminoGC_low/orthogroups.txt",
        "../../results/09_species_phylogeny_aminoGC/aminoGC_high/orthogroups.txt",
        "../../results/09_species_phylogeny_aminoGC/aminoGC.tsv",
        expand("../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.trimal.treefile", data=["aminoGC_high", "aminoGC_low"]),
        expand("../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.treefile", data=["aminoGC_high", "aminoGC_low"]),


rule alignment_aminogc:
    """
    To run this step download alignment_pruner.pl from
    https://github.com/novigit/davinciCode/blob/master/perl/alignment_pruner.pl
    and place it in the bin folder in the root of this directory.
    """
    input:
        "../../results/06_align_sco/alignments/{orthogroup}.aln"
    output:
        "../../results/09_species_phylogeny_aminoGC/alignment_pruner/alignments/{orthogroup}.tsv"
    conda:
        "../envs/alignment_prune.yaml"
    shell:
        "perl ../../bin/alignment_pruner.pl --file {input} --aminogc > {output}"

rule parse_results:
    input:
        expand("../../results/09_species_phylogeny_aminoGC/alignment_pruner/alignments/{orthogroup}.tsv", orthogroup=orthogroups)
    output:
        high_bias="../../results/09_species_phylogeny_aminoGC/aminoGC_high/orthogroups.txt",
        low_bias="../../results/09_species_phylogeny_aminoGC/aminoGC_low/orthogroups.txt",
        all_values="../../results/09_species_phylogeny_aminoGC/aminoGC.tsv"
    shell:
        "python parse_aminogc.py {input} {output.low_bias} {output.high_bias} {output.all_values}"

rule concatenate_aminoGC:
    input:
        protein_mappings="../../results/01_download_data/protein_mapping.tsv",
        alignment_dir="../../results/06_align_sco/trimmed_alignments",
        orthogroups="../../results/09_species_phylogeny_aminoGC/{data}/orthogroups.txt"
    output:
        concat="../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.trimal.aln"
    shell:
        "python ../concat_alignment.py --protein_mapping {input.protein_mappings} --orthogroups {input.orthogroups} --alignments {input.alignment_dir} --concatenated {output.concat}"

rule concatenate_aminoGC_untrimmed:
    input:
        protein_mappings="../../results/01_download_data/protein_mapping.tsv",
        alignment_dir="../../results/06_align_sco/alignments",
        orthogroups="../../results/09_species_phylogeny_aminoGC/{data}/orthogroups.txt"
    output:
        concat="../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.aln"
    shell:
        "python ../concat_alignment.py --protein_mapping {input.protein_mappings} --orthogroups {input.orthogroups} --alignments {input.alignment_dir} --concatenated {output.concat}"

rule species_phylogeny:
    input:
        concat="../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.trimal.aln",
    output:
        "../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.trimal.treefile",
    params:
        pre="../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.trimal",
        uf=1000,
        alrt=1000
    conda:
        "../envs/iqtree.yaml"
    threads:
        8
    shell:
        "iqtree -s {input.concat} -m LG -pre {params.pre} -bnni -bb {params.uf} -nt {threads} -redo"

rule species_phylogeny_untrimmed:
    input:
        concat="../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.aln",
    output:
        "../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.treefile",
    params:
        pre="../../results/09_species_phylogeny_aminoGC/{data}/trimal/{data}_orthogroups.",
        uf=1000,
        alrt=1000
    conda:
        "../envs/iqtree.yaml"
    threads:
        8
    shell:
        "iqtree -s {input.concat} -m LG -pre {params.pre} -bnni -bb {params.uf} -nt {threads} -redo"
