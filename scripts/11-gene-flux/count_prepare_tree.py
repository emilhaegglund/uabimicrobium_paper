#!/usr/bin/env python

import ete3
import sys
import re
import os
import argparse
import pandas as pd

def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', required=True)
    parser.add_argument('--outgroup', required=True)
    parser.add_argument('--tree_out', required=True)

    return parser.parse_args()

args = parse_command_line()

# Read tree
tree = ete3.Tree(args.tree, format=1)

# Search the fasta file name for the accession id for
# the outgroup genomes.
genomes_df = pd.read_csv(args.outgroup)
outgroup_genomes_df = genomes_df[genomes_df['Group'] == 'Outgroup']
outgroup = outgroup_genomes_df['Assembly Accession'].to_list()

for n in tree.traverse():

    if n.is_leaf():
        x = re.search('GCA_[0-9]*\.[0-9]', n.name)
        accession = x.group()

        if accession not in outgroup:
            first_root = n.name

tree.set_outgroup(first_root)
ancestor = tree.get_common_ancestor(outgroup)
tree.set_outgroup(ancestor)

print(tree)

# Rename nodes

node_counter = 0
for node in tree.traverse():
    node_name = 'N' + str(node_counter)

    if not node.is_leaf():  # Add names to internal nodes
        node.name = node_name
        print(node.name)

    else:  # Remove underscore from taxa name
        new_node_name = str(node.name).replace('_', '')
        node.name = new_node_name

    node_counter += 1

# Save to file

tree.write(format=1, outfile=args.tree_out)
