#!/usr/bin/env python
# from emils script

from ete3 import Tree, TextFace, TreeStyle
import pandas as pd
import argparse


def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    #parser.add_argument("--pkinase_og_result", required=True)
    parser.add_argument("--pkinase_species_result", required=True)
    parser.add_argument("--gained_og", required=True)
    parser.add_argument("--lost_og", required=True)
    parser.add_argument("--genome_stat", required=True)
    return parser.parse_args()


def read_tree(tree_file):
    """
    read the phylogeny
    """
    tree = Tree(tree_file, format=1)
    ts = TreeStyle()
    ts.show_leaf_name = False
    return tree, ts

def get_species(name, genome_stat):

    species = genome_stat[genome_stat['assembly_accession']== 'GCA_' + str(name.split('GCA')[-1])]['organism_name'].values.tolist()[0]
    return species

def get_group(name, genome_stat):
    group = genome_stat[genome_stat['assembly_accession']== 'GCA_' + str(name.split('GCA')[-1])]['group'].values.tolist()[0]
    return group

def render_tree(tree, gains_df, loss_df, pkinase_og, pkinase_species, genome_stat):
    """
    create a tree figure with the abscence and presence of pkinase domain found bu hmmer searc
    """
    node_counter = 0  # Keep track of the number of nodes visited for node naming purpose

    for n in tree.traverse():
        color = {'Isosphaeraceae':'#5FB3B3',
               'Pirellula':'#C0C5CE',
               'Gemmataceae':'#AB7968',
               'Saltatorellus':'#EC5F67',
               'Brocadiaceae':'#6699CC',
               'Bythopirellula':'#F7C863',
               'Gimesia':'#99C795',
               'Orphan':'#000000',
               'Outgroup':'#000000',
               'Uab':'#C594C5',
               'Phycisphaerae':'#F59157'}

        if not n.is_leaf():
            n.name = "N" + str(node_counter)  # Name the nodes in the same way as in the tree used in Count
        if n.name != "N0":

            ####### Need to fix so we read the parsimony files from count 
            gained_og = gains_df[gains_df['node']== n.name]['orthogroup'].values.tolist()
            lost_og = loss_df[loss_df['node']== n.name]['orthogroup'].values.tolist()

            ### extract the number of gained and lost pkinase ogs
            pkinase_og_gained = []
            for og in pkinase_og:
                if og in gained_og and og not in pkinase_og_gained:
                    pkinase_og_gained.append(og)

            pkinase_og_lost = []
            for og in pkinase_og:
                if og in lost_og and og not in pkinase_og_lost:
                    pkinase_og_lost.append(og)

            # Add gain and loss information to plot of species phylogeny

            n.add_face(
                TextFace(str(len(pkinase_og_gained)), fgcolor="green"),
                column=0,
                position="branch-top",
            )
            n.add_face(
                TextFace("/" + str(len(pkinase_og_lost)), fgcolor="red"),
                column=1,
                position="branch-top",
            )

            # Add node name to each node in the plot of the species phylogeny
            if n.is_leaf():
                species_name = get_species(n.name, genome_stat)
                group_name = get_group(n.name, genome_stat)
                if species_name == 'Planctomycetes bacterium SRT547':
                    species_name = 'Candidatus Uabimicrobium amorphum'
                    group_name='Uab'
                    print(pkinase_og_gained)

                n.add_face(TextFace(str(species_name), fgcolor=color[group_name]), column=0, position="branch-right")
    
                acc_in_node = 'GCA_' + str(n.name.split('GCA')[-1])
                #print(pkinase_species.columns)
                singletons = pkinase_species[pkinase_species['accession']==acc_in_node]['number of orthomcl singletons']
                all_pkinase = pkinase_species[pkinase_species['accession']==acc_in_node]['nr of pkinases']
                frac = pkinase_species[pkinase_species['accession']==acc_in_node]['fraction of pkinases (%)']

                if singletons.empty==False:
                    singletons = singletons.values.tolist()[0]
                else:
                    singletons = 0

                if all_pkinase.empty==False:
                    all_pkinase = all_pkinase.values.tolist()[0]
                else:
                    all_pkinase = 0

                if frac.empty==False:
                    frac = frac.values.tolist()[0]
                else:
                    frac = 0

                n.add_face(
                TextFace(str(singletons), fgcolor="black"), # blue
                column=0,
                position="branch-bottom",
                )

                n.add_face(
                TextFace("/" + str(all_pkinase), fgcolor="black"), #purple
                column=1,
                position="branch-bottom",
                )

                n.add_face(
                TextFace("/" + str(frac), fgcolor="black"),
                column=2,
                position="branch-bottom",
                )
            #else:
            #    n.add_face(TextFace(str(n.name)), column=0, position="branch-bottom")
        node_counter += 1
    return tree

def save_tree(tree, ts, out_file_name):
    """
    Plot species tree
    """
    f = tree.render(out_file_name, tree_style=ts)
    return f


def root_tree(tree, genome_stat):

    outgroup_genomes = genome_stat[genome_stat['group']=='Outgroup']['assembly_accession'].drop_duplicates().values.tolist()
    outgroup_genomes_mod = []

    for acc in outgroup_genomes:
        outgroup_genomes_mod.append(acc.replace('_', ''))

    common = tree.get_common_ancestor(outgroup_genomes_mod)
    tree.set_outgroup(common)
    #Root = tree.get_midpoint_outgroup()
    #tree.set_outgroup(Root)
    tree.ladderize(direction=1)
    
    for node in tree.traverse():
        if node.name == 'N1' or node.name == 'N2':
            node.ladderize(direction=0)
    return tree

def main(tree, pkinase_og_result, pkinase_species_result, gained_og, lost_og, genome_stat, out_file_name):
    
    ## extract input settings (paths and options)
    ##args = parse_command_line()

    # read files
    df_pkinase = pd.read_csv(pkinase_species_result, sep="\t")
    gained_og = pd.read_csv(gained_og, sep="\t")
    lost_og = pd.read_csv(lost_og, sep="\t")
    genome_stat = pd.read_csv(genome_stat, sep="\t")
    tree, ts = read_tree(tree)
    
    pkinase_og_list = df_pkinase['orthomcl OGs.1'].values.tolist()
    pkinase_og = []
    for l_sub in pkinase_og_list:
        l_sub= str(l_sub)
        if l_sub not in pkinase_og and l_sub != [] and l_sub!= 'nan':
            pkinase_og += l_sub.split(';')
        
    # create tree figure with pkinase presence and absence
    tree = render_tree(tree, gained_og, lost_og, pkinase_og, df_pkinase, genome_stat)
    tree = root_tree(tree, genome_stat)

    f = save_tree(tree, ts, out_file_name)

    return f

## input
tree = snakemake.input.tree
pkinase_og_result = '' #snakemake.input.pkinase_og_result
pkinase_species_result = snakemake.input.pkinase_species_result
gained_og = snakemake.input.gained_og
lost_og = snakemake.input.lost_og
genome_stat = snakemake.input.genome_stat
out_file_name = snakemake.output.out_file_name

# trun the script
main(tree, pkinase_og_result, pkinase_species_result, gained_og, lost_og, genome_stat, out_file_name)

