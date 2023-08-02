#!/usr/bin/env python
"""
Extract the orthogroups that have been gained and lost for each node in a phylogeny from
the output of Count parsimony method. The numbers for each node are plotted on the
phylogeny using ete. Also, two .tsv files are written, one for gains and one for losses,
containing annotation information for each protein together with which node they were
gained or lost in.
"""
from ete3 import Tree, TextFace, TreeStyle
import pandas as pd
import argparse
import os
import sys

def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', required=True)
    parser.add_argument('--count_result', required=True)
    parser.add_argument('--orthogroups', required=True)
    parser.add_argument('--clustering', required=True)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--outdir', required=True)

    return parser.parse_args()

args = parse_command_line()

# Read phylogeny
tree = Tree(args.tree, format=1)
families_df = pd.read_csv(args.count_result, sep='\t')

ts = TreeStyle()
ts.show_leaf_name = False

node_counter = 0  # Keep track of the number of nodes visited for node naming purpose

# Read orthogroups
orthogroups_dict = {}
with open(args.orthogroups, 'r') as f:
    for line in f:
        line = line.strip('\n')
        if args.clustering == 'orthofinder':
            line = line.split()
            orthogroup = line[0][:-1]
        elif args.clustering == "orthomcl":
            line = line.split()
            orthogroup = line[0]
        orthogroups = line[1:]
        orthogroups_dict[orthogroup] = orthogroups

# Dicts to store annotation data for each orthogroup
gains_dict = {
    "node": [],
    "orthogroup": [],
    "og_counter": [],
    "protein_id": [],
}

loss_dict = {
    "node": [],
    "orthogroup": [],
    "og_counter": [],
    "protein_id": [],
}
total_gains = 0
total_loss = 0
for n in tree.traverse():
    og_counter = 0
    if not n.is_leaf():
        n.name = 'N' + str(node_counter)
#    if verbose:
#        print(n.name)
#        print(n.up.name)
    if n.up != None:
        if n.up.name != 'N0':
#            print(n.name)
            diff = families_df[families_df[n.up.name] == 0][n.name] - families_df[families_df[n.up.name] == 0][n.up.name]
            sub_df = families_df[families_df[n.up.name] == 0]
            gains_df = sub_df[diff > 0][['name', n.name]].sort_values(by=n.name)
            gains_df['node'] = n.name

            for orthogroup in gains_df['name'].tolist():
                og_counter += 1
                for protein_id in orthogroups_dict[orthogroup]:
                    gains_dict['node'].append(n.name)
                    gains_dict['orthogroup'].append(orthogroup)
                    gains_dict['og_counter'].append(og_counter)
                    gains_dict['protein_id'].append(protein_id)

            n_family_gains = gains_df.shape[0]
            total_gains += n_family_gains

            n.add_face(TextFace(str(n_family_gains), fgcolor="green"), column=0, position='branch-top')

            diff = families_df[families_df[n.name] == 0][n.up.name] - families_df[families_df[n.name] == 0][n.name]
            sub_df = families_df[families_df[n.name] == 0]
            loss_df = sub_df[diff > 0][['name', n.name]].sort_values(by=n.name)
            loss_df['node'] = n.name

            for orthogroup in loss_df['name'].tolist():
                for protein_id in orthogroups_dict[orthogroup]:
                    loss_dict['node'].append(n.name)
                    loss_dict['orthogroup'].append(orthogroup)
                    loss_dict['protein_id'].append(protein_id)
                    loss_dict['og_counter'].append(og_counter)


            n_family_loss = loss_df.shape[0]
            total_loss += n_family_loss
            n.add_face(TextFace('/' + str(n_family_loss), fgcolor="red"), column=1, position='branch-top')

        if n.is_leaf():
            n.add_face(TextFace(str(n.name)), column=0, position='branch-right')
        else:
            n.add_face(TextFace(str(n.name)), column=0, position='branch-bottom')
    node_counter += 1

tree.render(args.prefix + '_gain_loss.tree.png', tree_style=ts)

all_gains_df = pd.DataFrame.from_dict(gains_dict)
all_gains_df[['node', 'orthogroup', 'og_counter', 'protein_id']].to_csv(args.prefix + '_gains.tsv', sep='\t', index=False)
#print(gains_dict)
all_loss_df = pd.DataFrame.from_dict(loss_dict)
all_loss_df[['node', 'orthogroup', 'og_counter', 'protein_id']].to_csv(args.prefix + '_loss.tsv', sep='\t', index=False)

