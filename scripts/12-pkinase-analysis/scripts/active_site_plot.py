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

def open_file(my_file):
    """
    open the hmm file as a dataframe
    """
    #col_names = ['target name','accession','tlen + query name','accession domain','qlen','E-value', 'score', 'bias','#', 'of', 'c-Evalue', 'i-Evalue', 'score2', 'bias2', 'domain from', 'domain to', 'protein from', 'protein to', 'from', 'to/acc/description of target']
    
    col_names = ['target name','accession','tlen', 'query name','accession domain','qlen','E-value', 'score', 'bias','#', 'of', 'c-Evalue', 'i-Evalue', 'score2', 'bias2', 'domain from', 'domain to', 'protein from', 'protein to']
    
    df= pd.read_csv(my_file, sep='\s',comment='#', header=None, names= col_names, skipinitialspace=True, usecols=col_names)
    
    df = df.fillna('Missing')

    return df

def open_TM(TM_file):
    """
    open the Transmembrane prediction file
    """
    TM_data = pd.read_csv(TM_file, header = None, sep='\t')
    return TM_data

############### Extract info from files ############################################



def extract_protein_ids(protein_data):
    """
    extract the protein ids
    """
    protein_ids = protein_data['target name']
    protein_ids = protein_ids.drop_duplicates()
    return protein_ids

def protein_length(protein_id, protein_data):
    """
    extract the protein length
    """

    p_length = protein_data[protein_data['target name'] == str(protein_id).rstrip()]['tlen'].drop_duplicates().values.tolist()[0]
    p_length = int(p_length)
    return p_length


def extract_domains(protein_id, protein_data):
    """
    extract the domains with domain id, start and stop for each domain prediction per protein
    """
    p_info = protein_data[protein_data['target name'] == str(protein_id).rstrip()]
    domain_info = p_info[['accession domain','protein from', 'protein to', 'domain from', 'domain to']]
    return domain_info


############### Figure text ############################################


def add_text_2_fig(grid, section, protein_name, p_length,space):
    """
    Add protein length and protein name as text to the images
    """

    plt.text(grid[section][0] + p_length + 100, grid[section][1],  str(p_length), fontsize = 10, horizontalalignment = 'center')   
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
    Add domains to the figure as orange rectangles
    """

    domain = str(domain_info[0]) # name
    index = [int(domain_info[1]), int(domain_info[2])] # start and stop
    rect = mpatches.Rectangle([grid[section][0]+index[0], grid[section][1]], (index[1] - index[0]), 50) #50 is the hight of the domain
    label([grid[section][0]+index[0], grid[section][1]-44], str(domain_info[3]) + '-' + str(domain_info[4]))
    patches.append(rect)
    

    colors.append('#8762aa')#[245/255.0,145/255.0,64/255.0]) # make the domains scilifelab-orange

    return patches, colors


############## Save figure with final layout settings ##################

def save_fig(protein_ids, out_name, patches, colors, plt, ax, x_lim, y_lim):
    """
    Save the figure with layout options
    """
    collection = PatchCollection(patches)
    collection.set_color(colors)
    ax.add_collection(collection)
    ax.figure.set_size_inches(60,0.5*int(len(protein_ids)))# 25 for pvc
    ax.set_xlim([-200,x_lim]) #1500
    ax.set_ylim([-2,y_lim]) #30000
    #plt.axis('off')
    ax.axis('off')
    ax.set_title(out_name.split('/')[-1].split('.')[0])
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    
    #plt.figure.frameon(False)
    plt.savefig(out_name, orientation = 'portrait', bbox_inches='tight')

    return plt



def create_protein_domain(sec, patches, colors, section, protein_id, interpro_df, grid):#, TM_data):
    """
    create a protein domain figure for one protein
    """
    domain_info = extract_domains(protein_id, interpro_df).values.tolist()

    p_length = protein_length(protein_id,interpro_df)
    p_amount = 0
    
    ## add protein id and protein length as text in the figure
    add_text_2_fig(grid, section, str(protein_id), p_length, 90)
    
    ## add a black line representing the protein length
    patches, colors = create_line(grid, section, p_length, patches, colors)

    ## add domain region in the protein

    for domain in domain_info:
        patches, colors = add_domains_2_fig(grid, section, domain, patches, colors)
        
    
    return patches, colors, section, grid


def generate_domain_fig(sec, interpro_df, protein_ids):#, TM_data):
    """
    Generate the domain figures for all proteins
    """
    ## figure variables
    fig, ax= plt.subplots()
    #ax = plt.axes()
    patches=[]
    colors=[]
    grid = np.mgrid[0:1,0:int(201 * len(protein_ids))].reshape(2,-1).T 
    #grid = np.mgrid[0:1:1j,0:4500:540j].reshape(2,-1).T # create grid to build the images on
    section = 0
                  

    # create domain figure for each protein in the interproscan file
    for protein_id in protein_ids:
        section+=sec
        
        patches, colors, section, grid = create_protein_domain(sec,patches, colors, section, protein_id, interpro_df, grid)#, TM_data)
    
                
    return patches, colors, fig, ax

def main(out_name, hmm, pid_list):
    """
    Generates a domain architecture figure of proteins in an interproscan file. 
    Use the input flag analysis to change between domain prediction softwares, default is Pfam.
    """

    hmm_df = open_file(hmm)
    
    if pid_list != []:
        hmm_df = hmm_df[hmm_df['target name'].isin(pid_list)]
    
    protein_ids = extract_protein_ids(hmm_df)
    l_max = 0

    for l in hmm_df['tlen'].drop_duplicates().values.tolist():
        if l > l_max:
            l_max=l 
    
    x_lim = l_max + 1500 #110
    y_lim = len(protein_ids)* 201 #240
    sec = 200
    patches, colors, plt, ax = generate_domain_fig(sec, hmm_df, protein_ids)

    f = save_fig(protein_ids, out_name, patches, colors, plt, ax, x_lim, y_lim)

    return f

############## Run script ##################
parser.add_argument("--out_file")
parser.add_argument("--hmm")
parser.add_argument("--pid_list")

args = parser.parse_args()


hmm = args.hmm
out_name = args.out_file
if args.pid_list:
    pid_list = pd.read_csv(str(args.pid_list), sep='\t')['protein_id'].values.tolist()
else:
    pid_list = []

main(out_name, hmm, pid_list)


