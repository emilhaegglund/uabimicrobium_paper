#! /usr/bin/env/python
"""
Outputs a concatenated alignment and a partition file.
"""
import os
from Bio import SeqIO
import re
import sys
import pandas as pd
import argparse


def read_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignments", required=True)
    parser.add_argument("--protein_mapping", "-m", required=True)
    parser.add_argument("--concatenated", required=True)
    parser.add_argument("--partition")
    parser.add_argument("--orthogroups")
    parser.add_argument("--species")
    parser.add_argument("--evo_models")

    return parser.parse_args()


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


def get_orthogroups_from_dir(alignment_dir):
    """
    Use all orthogroups that are located in the alignment directory.

    Return list of all orthogroups.
    """
    orthogroups = []
    for f in os.listdir(alignment_dir):
        orthogroup = f.split(".")[0]
        orthogroups.append(orthogroup)

    return orthogroups


def get_orthogroups_from_file(orthogroup_file_path):
    """
    Read orthogroups that will be included in the final alignment
    from a file.

    Return list of selected orthogroups.
    """
    orthogroups = []
    with open(orthogroup_file_path, "r") as f:
        for line in f:
            orthogroup = line.strip()
            orthogroups.append(orthogroup)

    return orthogroups


def get_species_from_file(species_file_path):
    """
    Read species that will be included in final alignment from file.

    Return a dictionary where each specie is a key and the value is
    placeholder for the alignment.
    """
    alignments = {}  # store all alignments
    with open(args.species, "r") as f:
        for line in f:
            specie = line.strip()
            all_species.append(specie)
            if specie not in alignments.keys():
                alignments[specie] = ""


def get_species_from_alignments(alignment_dir, protein_mappings, orthogroups):
    """
    Loop over all alignments to find all species that will be
    included in final alignment.

    Return a dictionary where each specie is a key and the value is
    placeholder for the alignment.
    """
    alignments = {}  # Store all alignments
    for f in os.listdir(alignment_dir):
        if f.split(".")[0] in orthogroups:
            for aln_record in SeqIO.parse(os.path.join(alignment_dir, f), "fasta"):
                specie = protein_mappings[aln_record.id]
                if specie not in alignments.keys():
                    alignments[specie] = ""  # Create a placeholder for the alignment

    return alignments


def concatenate_alignment(
    alignments,
    alignment_dir,
    protein_mappings,
    orthogroups,
    args
):
    """
    """
    if args.partition:
        partition = open(args.partition, "w")
    if args.evo_models:
        evo_models = pd.read_csv(args.evo_models, sep="\t")
        evo_models.rename(columns={"Unnamed: 0": "orthogroup"}, inplace=True)
    start = 0  # Start of the partition
    partition_counter = 1
    all_species = list(alignments.keys())
    # Concatenate all alignments
    for f in os.listdir(alignment_dir):
        orthogroup = f.split(".")[0]
        if orthogroup in orthogroups:
            aln_record = SeqIO.parse(os.path.join(alignment_dir, f), "fasta")
            covered_species = []
            for aln in aln_record:
                aln_length = len(aln.seq)
                specie = protein_mappings[aln.id]
                alignments[specie] += str(aln.seq)  # then concatenate the alignment
                covered_species.append(specie)
            missing_species = list(set(all_species) - set(covered_species))
            # insert deletions for missing genes
            for specie in missing_species:
                alignments[specie] += "-" * aln_length

            # Extract the evolutionary model to use for this orthogroup
            if args.evo_models:
                model = evo_models[evo_models["orthogroup"] == orthogroup][
                    ["subst_model"]
                ].values[0][
                    0
                ]  # evolution model to use
            else:
                model = "AUTO"

            if args.partition:
                partition.write(
                    model
                    + ", part"
                    + str(partition_counter)
                    + " = "
                    + str(start + 1)
                    + "-"
                    + str(start + aln_length)
                    + "\n"
                )
                start += aln_length  # Update the start position for the next partition
                partition_counter += 1

    if args.partition:
        partition.close()

    return alignments



# Read command line arguments and set the paths to use
args = read_command_line()
alignment_dir = args.alignments
protein_mapping_path = args.protein_mapping
protein_mappings = read_protein_mapping(protein_mapping_path)
if args.evo_models:
    print(evo_models)
concat_aln_path = args.concatenated


# Set which orthogroups to include in alignment
if args.orthogroups:
    orthogroups = get_orthogroups_from_file(args.orthogroups)
else:
    orthogroups = get_orthogroups_from_dir(alignment_dir)
print("Orthogroups included: " + str(len(orthogroups)))

# Set which speciees to include in alignment
if args.species:
    alignments = get_species_from_file(args.species)
if not args.species:
    alignments = get_species_from_alignments(
        alignment_dir, protein_mappings, orthogroups
    )
print("Species included: " + str(len(alignments.keys())))

alignments = concatenate_alignment(
    alignments, alignment_dir, protein_mappings, orthogroups, args
)

# Write the concatenated alignment to new file
concat_aln_file = open(concat_aln_path, "w")
for specie in alignments.keys():
    concat_aln_file.write(">" + specie + "\n")
    concat_aln_file.write(alignments[specie] + "\n")
concat_aln_file.close()
