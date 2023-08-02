rule all:
    input:
        "../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.orthomcl_COUNT_results.families.txt",

# 1. prepare the ML tree from figure 1, to be used for the gain-loss analysis. (The ML species tree from figure 1 is based on orthogroups created by the software orthoFinder)
rule count_prepare_tree:
     """
     Remove underscore from the taxa-names and add node-names
     """
     input:
         tree_in="../../results/07_species_phylogeny_single_copy_orthogroups/bmge/species_phylogeny.C60.alltrees_blength.contree.treefile",
         genomes_table="../../data/accession_numbers.csv"
     output:
         tree_out="../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.count.treefile"
     conda:
         "../envs/ete.yaml"
     shell:
         "python count_prepare_tree.py --tree {input.tree_in} --outgroup {input.genomes_table} --tree_out {output.tree_out}"

# 2.
rule count_prepare_table:
     """
     remove underscore from names
     """
     input:
         "../../results/03_orthomcl_clustering/orthomcl_results/Orthogroups.GeneCount.tsv"
     output:
         "../../results/11_gene_flux/bmge/orthomcl/03_orthomcl_Orthogroups.GeneCount.count.tsv"
     shell:
         "python count_prepare_table_orthomcl.py --protein_families {input} --out {output}"

rule count_parsimony:
    """
    To run this step, download Count from
    http://www.iro.umontreal.ca/~csuros/gene_content/count.html
    and place it in the bin-directory of this repository.
    """
    input:
         orthogroup_count="../../results/11_gene_flux/bmge/orthomcl/03_orthomcl_Orthogroups.GeneCount.count.tsv",
         tree="../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.count.treefile"
    output:
         "../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.orthomcl_COUNT_results.txt"
    params:
         gain=2
    shell:
         """
         java -Xmx16G -cp ../../bin/Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner -gain {params.gain} {input.tree} {input.orthogroup_count} > {output}
         """

rule get_family_changes:
     input:
         "../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.orthomcl_COUNT_results.txt"
     output:
         "../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.orthomcl_COUNT_results.families.txt"
     shell:
         "grep '# FAMILY' {input} > {output}"

rule get_gain_loss:
    input:
        tree = '../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.count.treefile',
        count_result = '../../results/11_gene_flux/bmge/orthomcl/08_species_phylogeny.C60.orthomcl_COUNT_results.families.txt',
        og = "../../results/03_orthomcl_clustering/orthomcl_results/Orthogroups.txt"
    output:
        "results/11_gene_flux/bmge/orthomcl_parsimony_gain_loss.tree.png",
        'results/11_gene_flux/bmge/orthomcl_parsimony_gains.tsv',
        'results/11_gene_flux/bmge/orthomcl_parsimony_loss.tsv'
    conda:
        "../envs/ete.yaml"
    shell:
        "python scripts/extract_gain_loss_parsimony.py --tree {input.tree} --count_result {input.count_result} --clustering 'orthomcl' --prefix 'results/orthomcl_parsimony' --orthogroups {input.og}"
