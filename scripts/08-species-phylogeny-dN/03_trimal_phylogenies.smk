dataset = ['dN_high', 'dN_low']

rule all:
    input:
        expand("../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.trimal.treefile", data=dataset),
        expand("../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.treefile", data=dataset)

rule concatenate_dN:
    input:
        protein_mappings="../../results/01_download_data/protein_mapping.tsv",
        alignment_dir="../../results/06_align_sco/trimmed_alignments",
        orthogroups="../../results/08_species_phylogeny_slow_evolving/{data}/orthogroups.txt"
    output:
        concat="../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.trimal.aln"
    shell:
        "python ../concat_alignment.py --protein_mapping {input.protein_mappings} --orthogroups {input.orthogroups} --alignments {input.alignment_dir} --concatenated {output.concat}"

rule concatenate_dN_untrimmed:
    input:
        protein_mappings="../../results/01_download_data/protein_mapping.tsv",
        alignment_dir="../../results/06_align_sco/alignments",
        orthogroups="../../results/08_species_phylogeny_slow_evolving/{data}/orthogroups.txt"
    output:
        concat="../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.aln"
    shell:
        "python ../concat_alignment.py --protein_mapping {input.protein_mappings} --orthogroups {input.orthogroups} --alignments {input.alignment_dir} --concatenated {output.concat}"

rule species_phylogeny_trimmed:
    input:
        concat="../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.trimal.aln",
    output:
        "../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.trimal.treefile",
    params:
        pre="../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.trimal",
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
        concat="../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.aln",
    output:
        "../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups.treefile",
    params:
        pre="../../results/08_species_phylogeny_slow_evolving/{data}/trimal/{data}_orthogroups",
        uf=1000,
        alrt=1000
    conda:
        "../envs/iqtree.yaml"
    threads:
        8
    shell:
        "iqtree -s {input.concat} -m LG -pre {params.pre} -bnni -bb {params.uf} -nt {threads} -redo"
