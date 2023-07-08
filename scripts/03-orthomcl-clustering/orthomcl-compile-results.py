#!/usr/bin/env python
"""
Script to convert orthoMCL output to be in the same form as the output
of OrthoFinder.
"""
import argparse
import os
import pandas as pd
from Bio import SeqIO
import sys
import re


def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument("--orthomcl", required=True)
    parser.add_argument("--species_dictionary", required=True)
    parser.add_argument("--proteome_dictionary", required=True)
    parser.add_argument("--protein_mapping", required=True)
    parser.add_argument("--proteome", required=True)
    parser.add_argument("--fraction", required=True, type=float)
    parser.add_argument("--output", default="orthomcl_results")

    return parser.parse_args()


def read_dictionary(dict_path):
    """
    Read the dictionary file that maps the four letter abbrevation used by
    orthoMCL to the full species name.
    """
    species_dict = {}
    with open(dict_path, "r") as dict_file:
        for line in dict_file:
            line = line.strip("\n")
            line = line.split("\t")
            species_dict[line[0]] = line[1]

    return species_dict


def create_orthogroups(mcl_path):
    """
    Create a dictionary with all Orthogroups
    """
    orthogroup_dict = {}  # Dictionary to store orthogroups
    og_counter = 0
    with open(mcl_path, "r") as mcl_file:
        n = len(mcl_file.readlines())

    n = len(list(str(n))) + 2
    with open(mcl_path, "r") as mcl_file:
        for line in mcl_file:
            og_name = "OG" + str(og_counter).zfill(n)
            orthogroup_dict[og_name] = []
            line = line.strip("\n")
            line = line.split(" ")
            for gene in line[1:]:
                orthogroup_dict[og_name].append(gene)
            og_counter += 1

    return orthogroup_dict


def read_protein_mapping(protein_mapping_path):
    """
    Read file containing the mapping of record ids to assembly
    accessions.

    Return dictionary with the protein ids as keys.
    """
    protein_mappings = {}
    with open(protein_mapping_path, "r") as f:
        for line in f:
            line = line.strip()
            line = line.split("\t")
            protein_mappings[line[0]] = line[1]

    return protein_mappings


def create_gene_counts(orthogroup_dict, species_dict, protein_mapping, dir_path):
    """
    Create the GeneCount file
    """
    output_path = os.path.join(dir_path, "Orthogroups.GeneCount.tsv")
    species_list = [species_dict[specie] for specie in species_dict.keys()]
    trimmed_names = []
    for species in species_list:
        print(species)
        trimmed_names.append(os.path.splitext(species)[0])

    orthogroup_count = {}
    for orthogroup in orthogroup_dict.keys():
        og = orthogroup_dict[orthogroup]
        orthogroup_count[orthogroup] = {}
        species_in_og = []
        for gene in og:
            species_in_og.append(protein_mapping[gene])
        for specie in trimmed_names:
            orthogroup_count[orthogroup][specie] = species_in_og.count(specie)

    # Create row indexes similar to the ones in OrthoFinder
    # i.e. OG0012233
    ogs = len(list(orthogroup_dict.keys()))
    n = len(list(str(ogs))) + 2
    index = []
    for i in range(ogs):
        index.append("OG" + str(i).zfill(n))

    # Convert orthogroup dictionary to DataFrame for saving
    orthogroup_gene_count_df = pd.DataFrame(orthogroup_count)
    orthogroup_gene_count_df = orthogroup_gene_count_df.T
    orthogroup_gene_count_df.rename(species_dict, axis="columns", inplace=True)
    orthogroup_gene_count_df.to_csv(output_path, sep="\t")

    return orthogroup_gene_count_df


def write_single_copy_genes(orthogroups_count_df, dir_path):
    output_path = os.path.join(dir_path, "SingleCopyOrthogroups.txt")
    output = open(output_path, "w")
    count_groups = 0
    for index, row in orthogroups_count_df.iterrows():
        if all(count == 1 for count in row.tolist()):
            count_groups += 1
            output.write(row.name + "\n")
    output.close()
    print("Found " + str(count_groups) + " single copy orthogroups")


def write_almost_single_copies(orthogroups_count_df, dir_path, species_dict, fraction):
    """
    Find orthogroups that are single copy in at least a set
    fraction of the species.
    """
    output_path = os.path.join(dir_path, "AlmostSingleCopyOrthogroups.txt")
    output = open(output_path, "w")
    count_groups = 0
    for index, row in orthogroups_count_df.iterrows():
        one_zero_row = [
            i for i in row.tolist() if i in [0, 1]
        ]  # Extract OG with 0 or 1 as copy number
        ones = one_zero_row.count(1)
        if (
            len(one_zero_row) == len(row) and ones / len(row) > fraction
        ):  # Check that all species are either 1 or 0
            count_groups += 1
            output.write(row.name + "\n")
    output.close()
    print(
        "Found "
        + str(count_groups)
        + " single copy orhtogroups covering "
        + str(fraction)
        + " of the taxa"
    )


def find_species_specific_orthogroups(orthogroups_count_df, dir_path):
    species_specific = 0
    nr_species_specific_genes = 0
    all_species = 0
    for index, row in orthogroups_count_df.iterrows():
        zeros = 0
        genes = 0
        row_length = len(row)
        for count in row.tolist():
            if count == 0:
                zeros += 1
            else:
                genes = count
        if zeros == row_length - 1:
            nr_species_specific_genes += genes
            species_specific += 1
        elif zeros == 0:
            all_species += 1

    return species_specific, nr_species_specific_genes, all_species


def translate_proteinid(orthogroup_dict, protein_dict):
    """
    Rename the proteins in the orthogroups to their
    original ids.
    """
    genes_in_orthogroup = 0
    for orthogroup in orthogroup_dict:
        for i, protein in enumerate(orthogroup_dict[orthogroup]):
            old_id = protein_dict[protein]
            orthogroup_dict[orthogroup][i] = old_id
            genes_in_orthogroup += 1

    return orthogroup_dict, genes_in_orthogroup


def create_singletons(orthogroup_dict, protein_dict):
    """
    Place unassigned genes orhtogroups
    """
    proteins_in_orthogroup = []
    n_orthogroups = 0
    for orthogroup, proteins in orthogroup_dict.items():
        proteins_in_orthogroup += proteins
        n_orthogroups += 1

    all_proteins = []
    for new_id, protein in protein_dict.items():
        all_proteins.append(protein)

    singletons = set(all_proteins) - set(proteins_in_orthogroup)

    singletons_dict = {}
    for protein in singletons:
        singletons_dict[n_orthogroups] = [protein]
        n_orthogroups += 1

    return singletons_dict


def write_orthogroups(orthogroups, path):
    output_path = os.path.join(path, "Orthogroups.txt")
    output = open(output_path, "w")
    for orthogroup in orthogroups.keys():
        new_line = orthogroup + " " + " ".join(orthogroups[orthogroup]) + "\n"
        output.write(new_line)
    output.close()


def create_statistics(orthogroups, number_of_genes, orthogroups_count_df, dir_path):
    output_path = os.path.join(dir_path, "Statistics_Overall.csv")
    with open(output_path, "w") as output:
        output.write("Number of genes " + str(number_of_genes) + "\n")
        output.write(
            "Number of genes in orthogroups " + str(genes_in_orthogroup) + "\n"
        )
        output.write(
            "Number of unassigned genes "
            + str(number_of_genes - genes_in_orthogroup)
            + "\n"
        )
        output.write("Number of orthogroups " + str(len(orthogroups.keys())) + "\n")
        (
            species_specific,
            nr_species_specific_genes,
            all_species,
        ) = find_species_specific_orthogroups(orthogroups_count_df, dir_path)
        output.write(
            "Number of species-specific orthogroups " + str(species_specific) + "\n"
        )
        output.write(
            "Number of genes in species-specific orthogroups "
            + str(nr_species_specific_genes)
            + "\n"
        )
        output.write(
            "Mean orthogroup size "
            + str(genes_in_orthogroup / len(orthogroups.keys()))
            + "\n"
        )
        output.write(
            "Number of orthogroups with all species present " + str(all_species) + "\n"
        )


def create_sequence_dir(orthogroups, protein_dict, proteome_path, output_path):
    """
    For each orthogroup, create a fasta file containing the protein
    sequences included in the orthogroup.
    """
    fasta_dict = {}
    fastas = [
        fasta
        for fasta in os.listdir(proteome_path)
        if os.path.isfile(os.path.join(proteome_path, fasta))
    ]
    for fasta in fastas:
        fasta_path = os.path.join(proteome_path, fasta)
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id not in fasta_dict.keys():
                original_name = protein_dict[record.id]
                record.id = original_name
                fasta_dict[original_name] = record
            else:
                sys.exit("Error")

    number_of_genes = len(fasta_dict.keys())

    # Create the same directory structure as in OrthoFinder
    orthologoues_path = os.path.join(output_path, "Orthologues")
    sequence_path = os.path.join(orthologoues_path, "Sequences")
    if not os.path.exists(orthologoues_path):
        os.mkdir(orthologoues_path)
        os.mkdir(sequence_path)
    else:
        if not os.path.exists(sequence_path):
            os.mkdir(sequence_path)

    for orthogroup in orthogroups.keys():
        fasta_records = []
        output = os.path.join(
            output_path, "Orthologues", "Sequences", orthogroup + ".faa"
        )
        for gene in orthogroups[orthogroup]:
            fasta_records.append(fasta_dict[gene])
        SeqIO.write(fasta_records, output, "fasta")

    return number_of_genes


if __name__ == "__main__":
    args = parse_command_line()
    mcl_path = args.orthomcl
    species_dict = read_dictionary(args.species_dictionary)
    protein_dict = read_dictionary(args.proteome_dictionary)
    protein_mapping = read_protein_mapping(args.protein_mapping)
    print("Create orthogroups")
    orthogroups = create_orthogroups(args.orthomcl)
    orthogroups, genes_in_orthogroup = translate_proteinid(orthogroups, protein_dict)
    print("Write Orthogroups.txt")
    write_orthogroups(orthogroups, args.output)
    print("Write orhtogroups.genecount")
    orthogroups_count_df = create_gene_counts(
        orthogroups, species_dict, protein_mapping, args.output
    )
    print("Write singlecopygenes")
    write_single_copy_genes(orthogroups_count_df, args.output)
    print("Write Almost single copies")
    write_almost_single_copies(
        orthogroups_count_df, args.output, species_dict, fraction=args.fraction
    )
    print("Write Orthogroups to fasta")
    number_of_genes = create_sequence_dir(
        orthogroups, protein_dict, args.proteome, args.output
    )
    create_statistics(orthogroups, number_of_genes, orthogroups_count_df, args.output)
