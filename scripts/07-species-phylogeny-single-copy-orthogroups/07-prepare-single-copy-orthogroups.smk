rule all:
    input:
        "../../results/07_species_phylogeny_single_copy_orthogroups/orthogroups.txt"

rule select_orthogroups:
    input:
        "../../results/02_orthofinder/single_copy_orthogroups.tsv"
    output:
        "../../results/07_species_phylogeny_single_copy_orthogroups/orthogroups.txt"
    shell:
        "awk -F'\\t' '{{ print $2 }}' {input} > {output}"
