# Run the Slow evolving protein workflow, generating a phylogenetic tree and statistics

# 1. Extract the corresponding genes to the proteins in the single copy orthogroups.
# 2. Create MSA files of the gene sequence orthogroup files
# 3. Execute codeML to calculate dN and dS pairwise values, see the codeml.ctl file for settings
# 4. Parse the output from codeml to a file format easier to work with
# 5. Statistical analysis, generate a statistical table, a dN mean figure and Orthogroups with lowest mean dN are saved to a file to use for the the tree in step 6
# 6. Create the slow evolving phylogeny


import pandas as pd

df = pd.read_csv('../../results/02_orthofinder/single_copy_orthogroups.tsv', sep='\t', names=['n_species', 'orthogroup'])
orthogroups = df['orthogroup'].to_list()

rule all:
    input:
        expand('../../results/08_species_phylogeny_slow_evolving/orthogroup_nt_sequences/{orthogroup}.fna', orthogroup=orthogroups),
        expand('../../results/08_species_phylogeny_slow_evolving/orthogroup_nt_alignments/{orthogroup}.aln', orthogroup=orthogroups),
        expand('../../results/08_species_phylogeny_slow_evolving/CODEML/aa_msa_converted_to_nt/{orthogroup}.pal2nal', orthogroup=orthogroups),
        expand('../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/{orthogroup}.txt', orthogroup=orthogroups),
        expand('../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/2ML.dN', orthogroup=orthogroups),
        expand('../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/2ML.dS', orthogroup=orthogroups),
        expand('../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_parsed/{orthogroup}_dnds_table.tsv', orthogroup=orthogroups),
        '../../results/08_species_phylogeny_slow_evolving/dN_high/orthogroups.txt',
        '../../results/08_species_phylogeny_slow_evolving/dN_low/orthogroups.txt',
        '../../results/08_species_phylogeny_slow_evolving/dN_orthogroups.png',
        '../../results/08_species_phylogeny_slow_evolving/dNdS_statistics.tsv'


rule extract_OG_nt:
    input:
        '../../results/06_align_sco/alignments/{orthogroup}.aln',
        '../../results/01_download_data/genbank_genes/',
        '../../results/01_download_data/protein_mapping.tsv'
    output:
        '../../results/08_species_phylogeny_slow_evolving/orthogroup_nt_sequences/{orthogroup}.fna'
    script:
        'extract_OG_nt_seq.py'

rule msa_OG_nt:
    input:
        '../../results/08_species_phylogeny_slow_evolving/orthogroup_nt_sequences/{orthogroup}.fna'
    output:
        '../../results/08_species_phylogeny_slow_evolving/orthogroup_nt_alignments/{orthogroup}.aln'
    conda:
        "../envs/mafft.yaml"
    threads:
        4
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"

rule pal2nal:
    input:
        '../../results/06_align_sco/alignments/{orthogroup}.aln',
        '../../results/08_species_phylogeny_slow_evolving/orthogroup_nt_alignments/{orthogroup}.aln'
    output:
        '../../results/08_species_phylogeny_slow_evolving/CODEML/aa_msa_converted_to_nt/{orthogroup}.pal2nal'
    conda:
        "../envs/pal2nal.yaml"
    shell:
        'pal2nal.pl {input[0]} {input[1]} -output paml -nogap > {output}'


rule run_codeml:
    input:
        '../../results/08_species_phylogeny_slow_evolving/CODEML/aa_msa_converted_to_nt/{orthogroup}.pal2nal'
    output:
        '../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/{orthogroup}.txt',
        '../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/2ML.dN',
        '../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/2ML.dS'
    params:
        outdir = '../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/'
    conda:
        "../envs/paml.yaml"
    script:
        'run_codeml_biopython.py'

rule run_parser:
    input:
        '../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/2ML.dN',
        '../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_raw/{orthogroup}/2ML.dS'
    output:
        '../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_parsed/{orthogroup}_dnds_table.tsv'
    script:
        'create_dnds_table.py'

rule statistics:
    input:
        expand('../../results/08_species_phylogeny_slow_evolving/CODEML/codeml_parsed/{orthogroup}_dnds_table.tsv', orthogroup=orthogroups)
    output:
        '../../results/08_species_phylogeny_slow_evolving/dN_low/orthogroups.txt',
        '../../results/08_species_phylogeny_slow_evolving/dN_high/orthogroups.txt',
        '../../results/08_species_phylogeny_slow_evolving/dN_orthogroups.png',
        '../../results/08_species_phylogeny_slow_evolving/dNdS_statistics.tsv'
    script:
        'dnds_analysis.py'
