import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.colors as c
import pandas as pd
import numpy as np
import argparse
parser = argparse.ArgumentParser()

############### Open files ############################################

def open_file(my_file, general_domains):
    """
    open the interproscan file as a dataframe
    """
    interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession',
                           'Signature description','Start location', 'Stop location', 'Score', 'Status',
                           'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
                           'optional: GO terms', 'optional: Pathway annotation']

    my_header = ['Architecture_id','Sequence length','Signature accession','Start location', 'Stop location', 'Abundance in interproscan file', 'Protein_id', 'Signature description']

    if general_domains == 'No':
        df = pd.read_csv(my_file, header = None, sep='\t', names = interproscan_header)
    else:
        df= pd.read_csv(my_file, header = None, sep='\t', names = my_header)

    df = df.fillna('Missing')
    return df

def open_TM(TM_file):
    """
    open the Transmembrane prediction file
    """
    TM_data = pd.read_csv(TM_file, header = None, sep='\t')
    return TM_data

############### Extract info from files ############################################

def extract_analysis(domain_file, analysis, general_domains):
    """
    extract the domains of the interproscan file, default the analysis is Pfam
    """

    interpro_df = open_file(domain_file, general_domains)

    if general_domains == 'No':
        domain_df = interpro_df[interpro_df['Analysis']==analysis]
    else:
        domain_df = interpro_df

    return domain_df


def extract_protein_ids(protein_data):
    """
    extract the protein ids
    """
    protein_ids = protein_data['Architecture_id']
    protein_ids = protein_ids.drop_duplicates()
    return protein_ids

def protein_length(protein_id, protein_data):
    """
    extract the protein length
    """
    p_length = protein_data[protein_data['Architecture_id'] == str(protein_id).rstrip()].drop_duplicates()['Sequence length']
    p_length = int(p_length.values.tolist()[0])
    return p_length


def protein_amount(protein_id, protein_data):
    """
    extract the protein length
    """
    amount_s = protein_data[protein_data['Architecture_id'] == str(protein_id).rstrip()].drop_duplicates()['Abundance in interproscan file'].values.tolist()[0]
    return amount_s

def extract_domains(protein_id, protein_data):
    """
    extract the domains with domain id, start and stop for each domain prediction per protein
    """
    p_info = protein_data[protein_data['Architecture_id'] == str(protein_id).rstrip()]
    domain_info = p_info[['Signature accession','Start location', 'Stop location']]
    return domain_info

def extract_TM(protein_id, TM_data):
    """
    extract the transmembrane info / signaling too?
    """
    TM = TM_data[TM_data['protein_id'] == str(protein_id).rstrip()]
    return TM

def get_groups(interpro_df, protein_id, groups):
    group_list = []
    group_temp = ''

    real_protein_ids = interpro_df[interpro_df['Architecture_id']==protein_id]['Protein_id'].values.tolist()[0]

    for protein in real_protein_ids.split(';'):
        group_temp = groups[groups['pkinase protein ids'].str.contains(protein)]['group'].drop_duplicates().values.tolist()[0]

        if group_temp not in group_list:
            group_list.append(group_temp)
        group_temp = ''
    return group_list

############### Figure text ############################################


def add_text_2_fig(general_domain, grid, section, protein_name, p_length, p_amount, space):
    """
    Add protein length and protein name as text to the images
    """
    if general_domains == 'No':
        plt.text(grid[section][0] + p_length + 100, grid[section][1],  str(p_length), fontsize = 10, horizontalalignment = 'center')
    else:
        plt.text(grid[section][0] + p_length + 100, grid[section][1],  str(p_amount), fontsize = 10, horizontalalignment = 'center')
    plt.text(grid[section][0] - space, grid[section][1], protein_name + '  ', fontsize = 10, horizontalalignment = 'right')
    return plt

def label(xy,text):
    """
    Create the text label layout
    """
    plt.text(xy[0], xy[1] -15, text, ha='left', va='center',horizontalalignment= 'right',rotation =0, family='sans-serif', size = 8)
    return plt

def label_bold(xy,text):
    """
    Create the bold text label layout
    """
    plt.text(xy[0], xy[1] -15, text, ha='left', va='center', family='sans-serif', size = 8, weight='bold')
    return plt

def create_lable(grid, section, lable, patches, margin):
    """
    create lable text box of color and domain ids
    """
    plt.text(grid[section][0] + 2500, grid[section][1] - margin, str(lable), ha='left', va='center', family='sans-serif', size = 16) #, weight='bold')

    return plt

def create_lable_rect(grid, section, lable, color, patches, colors, margin):
    rect = mpatches.Rectangle([grid[section][0]+2400, grid[section][1]- margin - 25], 50, 50) #50 is the hight of the domain
    patches.append(rect)
    colors.append(color)
    return patches, colors

def create_lable_circle(grid, section, lable, color, patches, colors, margin):
    circle = mpatches.Ellipse((int(grid[section][0] + 2425), int(grid[section][1] - margin)), 25, 100)
    patches.append(circle)
    colors.append(color)
    return patches, colors
############# Add shapes; lines and rectangles representing protein length, domains and transmembrane regions

def create_line(grid, section, p_length, patches, colors):
    """
    create the line representing the protein length
    """
    rect = mpatches.Rectangle([grid[section][0],grid[section][1]+25],p_length,5) # create the black line that symbolises sequence length (grid[section][1]+25)
    patches.append(rect)
    colors.append([0.0,0.0,0.0])

    return patches, colors

def add_domains_2_fig(grid, section, domain_info, patches, colors):
    """
    Add domains to the figure as colored rectangles
    """
    top_10_domain_id_list_keys = top_10
    global colors_for_top_10_values
    colors_for_top_10_values = [  '#FFBF65', '#8DD78F', '#FF96C5','#C05780','#E77577', '#00B0BA', '#FFD872', '#4DD091', '#6C88C4', '#FF828B']
    color_dict = dict(zip(top_10_domain_id_list_keys, colors_for_top_10_values))

    domain = str(domain_info[0]) # name
    index = [int(domain_info[1]), int(domain_info[2])] # start and stop
    rect = mpatches.Rectangle([grid[section][0]+index[0], grid[section][1]], (index[1] - index[0]), 50) #50 is the hight of the domain
    label([grid[section][0]+index[0], grid[section][1]-44], domain)
    patches.append(rect)

    if domain not in ['SIGNAL', 'c_TM_n', 'n_TM_c']:

        if domain not in color_dict.keys():
            colors.append('#00A5E3')
        else:
            colors.append(color_dict[domain])

    else:
        colors.append('grey')

    return patches, colors

def add_group_circle_2_fig(grid, section, group, patches, colors, p_length, cnt_group):
    """
    add what group has this domain architecture
    """
    global color_dict
    color_dict = {'Isosphaeraceae':'#5FB3B3',
               'Pirellula':'#C0C5CE',
               'Gemmataceae':'#AB7968',
               'Saltatorellus':'#EC5F67',
               'Brocadiaceae':'#6699CC',
               'Bythopirellula':'#F7C863',
               'Gimesia':'#99C795',
               'Orphan':'#3f729b',
               'Outgroup':'#000000',
               'Orphan_Uab':'#C594C5',
               'Phycisphaerae':'#F59157',
               'Ascomycota':'#CC79A7',
               'Nematoda':'#3f729b',
                 'Proteobacteria':'#CC79A7'}

    buffer = 50*cnt_group
    circle = mpatches.Ellipse((int(grid[section][0]+ p_length  + 100  + int(buffer)), int(grid[section][1])),25,100)
    patches.append(circle)
    colors.append(color_dict[group])

    return patches, colors

def add_TM_2_fig(TM, section, grid, patches, colors):
    """
    Add transmembrane regions to the figure as blue rectangles
    """
    rect = mpatches.Rectangle([grid[section][0] + int(TM[1]), grid[section][1]], (int(TM[2])-int(TM[1])), 50)
    patches.append(rect)
    colors.append([17/255.0,134/255.0,195/255.0])

    return patches, colors

def add_signal_2_fig(SIG, section, grid, patches, colors):
    """
    Add transmembrane regions to the figure as blue rectangles
    """
    rect = mpatches.Rectangle([grid[section][0] + int(SIG[1]), grid[section][1]], (int(SIG[2])-int(SIG[1])), 50)
    patches.append(rect)
    colors.append([17/255.0,134/255.0,195/255.0]) # change color to green

    return patches, colors

############## Save figure with final layout settings ##################

def save_fig(protein_ids, out_name, patches, colors, plt, ax, x_lim, y_lim):
    """
    Save the figure with layout options
    """
    collection = PatchCollection(patches)
    collection.set_color(colors)
    ax.add_collection(collection)
    #size = 60
    #if out_name == 'results/pkinase_species_schematic_figure_other.pdf':
    #    size=60
    ax.figure.set_size_inches(90,0.5*int(len(protein_ids)+40))# 25 for pvc # 0.5 ### 60,0.7
    ax.set_xlim([-200,x_lim]) #1500
    ax.set_ylim([-2,y_lim]) #30000
    ax.axis('off')
    plt.savefig(out_name, orientation = 'portrait', bbox_inches='tight')

    return plt



def create_protein_domain(sec, general_domain, patches, colors, section, protein_id, interpro_df, grid, groups):#, TM_data):
    """
    create a protein domain figure for one protein
    """
    domain_info = extract_domains(protein_id, interpro_df).values.tolist()

    p_length = protein_length(protein_id, interpro_df)
    p_amount = 0

    if general_domain == 'Yes':
        p_amount = protein_amount(protein_id, interpro_df)

    ## add protein id and protein length as text in the figure
    add_text_2_fig(general_domain, grid, section, str(protein_id), p_length, p_amount, 90)

    ## add a black line representing the protein length
    patches, colors = create_line(grid, section, p_length, patches, colors)

    ## add domain region in the protein

    for domain in domain_info:
        if 'X' != domain[0]:
            patches, colors = add_domains_2_fig(grid, section, domain, patches, colors)

    ## add the group as a colored circle:
    cnt_group = 1
    for group in groups:
        patches, colors = add_group_circle_2_fig(grid, section, group, patches, colors, p_length, cnt_group)
        cnt_group+=1

    return patches, colors, section, grid

def get_anno(top_10, interpro_df):
    anno_list = []
    for domain_id in top_10:
        anno_list.append(interpro_df[interpro_df['Signature accession']==domain_id]['Signature description'].values.tolist()[0])
    return anno_list

def generate_domain_fig(sec, interpro_df, protein_ids, general_domain, groups, data_set):
    """
    Generate the domain figures for all proteins
    """
    ## figure variables
    fig, ax= plt.subplots()
    patches=[]
    colors=[]
    grid = np.mgrid[0:1,0:int(201 * len(protein_ids))].reshape(2,-1).T # create grid to build the images on
    #grid = np.mgrid[0:1:1j,0:4500:540j].reshape(2,-1).T # create grid to build the images on

    section = 0

    top_10_domain_colors = interpro_df
    top_10_domain_colors = top_10_domain_colors[top_10_domain_colors['Signature accession']!='Missing']
    top_10_domain_colors = top_10_domain_colors[['Signature accession']].value_counts().reset_index()
    top_10_domain_colors.columns = ['Signature accession', 'Count']
    top_10_domain_colors = top_10_domain_colors[~top_10_domain_colors['Signature accession'].isin(['SIGNAL', 'X', 'c_TM_n', 'n_TM_c'])]
    global top_10
    top_10 = top_10_domain_colors.iloc[0:10,0].values.tolist()
    top_10_anno = get_anno(top_10, interpro_df)


    # create domain figure for each protein in the interproscan file
    for protein_id in protein_ids:
        section+=sec

        # if schematic table as input:
        group_list = get_groups(interpro_df, protein_id, groups)

        patches, colors, section, grid = create_protein_domain(sec, general_domain, patches, colors, section, protein_id, interpro_df, grid, group_list)

    margin = 0
    counter = 0


    for lable in top_10_anno:
        create_lable(grid, section, str(lable), patches,margin)
        patches, colors = create_lable_rect(grid, section, lable, colors_for_top_10_values[counter], patches, colors, margin)
        margin += 150
        counter+=1

    # create a lable for the domains not in the top 10, that all has the same color:
    create_lable(grid, section, 'Other domains', patches,margin)
    patches, colors = create_lable_rect(grid, section, 'Other domains', '#00A5E3', patches, colors, margin)
    margin += 300

    if data_set == 'pvc':
        groups = groups[~groups['group'].isin(['Nematoda', 'Ascomycota'])]
    if data_set == 'other':
        groups = groups[groups['group'].isin(['Nematoda', 'Ascomycota'])]

    for group in groups['group'].drop_duplicates().values.tolist():
        create_lable(grid, section, str(group), patches, margin)
        patches, colors = create_lable_circle(grid, section, str(group), color_dict[group], patches, colors, margin)
        margin += 150

    return patches, colors, fig, ax

def main(Interproscan, analysis, out_name, general_domains, groups, data_set):
    """
    Generates a domain architecture figure of proteins in an interproscan file.
    Use the input flag analysis to change between domain prediction softwares, default is Pfam.
    """

    interpro_df = extract_analysis(Interproscan, analysis, general_domains)
    protein_ids = extract_protein_ids(interpro_df)

    x_lim = int(max(interpro_df['Sequence length'])) + 1500 + 500 #+ 500
    y_lim = len(protein_ids)* 201 #500#201#240
    sec = 200
    patches, colors, plt, ax = generate_domain_fig(sec, interpro_df, protein_ids, general_domains, groups, data_set)

    f = save_fig(protein_ids, out_name, patches, colors, plt, ax, x_lim, y_lim)

    return f

############## Run script ##################


global data_set

Interproscan = str(snakemake.input.schematic_table)
out_name = str(snakemake.output.out_file)
analysis = str(snakemake.params.analysis)
general_domains = 'Yes'
data_set = str(snakemake.wildcards.dataset)
groups_protein_ids = str(snakemake.input.groups_protein_ids)
groups_protein_ids = pd.read_csv(groups_protein_ids, sep='\t')[['group','pkinase protein ids']]

main(Interproscan, analysis, out_name, general_domains, groups_protein_ids, data_set)


