orthogroups = glob_wildcards('../../results/06_align_sco/alignments/{orthogroup}.aln').orthogroup

rule all:
    input:
        "../../results/09_species_phylogeny_aminoGC/aminoGC_low/orthogroups.txt",
        "../../results/09_species_phylogeny_aminoGC/aminoGC_high/orthogroups.txt",
        "../../results/09_species_phylogeny_aminoGC/aminoGC.tsv"

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
