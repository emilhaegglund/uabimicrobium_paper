import seaborn as sns
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib import rcParams
import argparse
import numpy as np
import matplotlib.colors

#def create_colormap():
#    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#FFFFE0","#9999CC"])
#    return cmap


def open_file(my_file):
    """
    open the interproscan file as a dataframe
    """
    interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession',
                           'Signature description','Start location', 'Stop location', 'Score', 'Status',
                           'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
                           'optional: GO terms', 'optional: Pathway annotation']
    df = pd.read_csv(my_file, header = None, sep='\t', names = interproscan_header)
    df = df.fillna('Missing')
    df = df[df['Analysis']=='Pfam']
    #pids = df[df['Signature accession'].isin(['PF00069'])]['protein_id'].drop_duplicates().values.tolist()
    #df = df[df['protein_id'].isin(pids)]
    df = df.sort_values(by=['protein_id','Start location'])

    return df


def count_abundance(domain_file):
    cnt_dom = pd.DataFrame(domain_file.value_counts()).reset_index()
    print(cnt_dom.columns)
    cnt_dom.columns = ['group','organism_name', 'Signature accession merged', 'Signature description merged', 'Count']
    #cnt_dom.to_csv('../results/heat_map_domain_abundance_dataset.tsv', sep='\t', index=False)
    cnt_dom = cnt_dom[['group','organism_name','Signature description merged', 'Count']]
    return cnt_dom

def absence_precence_table(cnt_dom):
    cnt_abs_pres = cnt_dom
    cnt_abs_pres['Count'].values[cnt_abs_pres['Count'] > 1.0] = 1.0
    cnt_abs_pres.drop_duplicates()
    return cnt_abs_pres

def cnt_dom_matrix(cnt_df):
    cnt_dom_matrix = cnt_df.pivot_table(index=['group', 'organism_name'], columns='Signature description merged',values='Count').fillna(0)
    cnt_dom_matrix = cnt_dom_matrix.T
    return cnt_dom_matrix


def create_heatmap(cnt_dom_matrix, outfile):

    sns.set(font_scale=1) #6699CC #C594C5
    c = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#FFFFE0","#6699CC"]) #"#FFFFE0","#9999CC" # "#3274A1","#E1812C"
    ax = sns.clustermap(cnt_dom_matrix, cmap=c, figsize=(20, 80), dendrogram_ratio=(.1, .5)) #"YlGnBu" 6.4, 4.8 # 20, 80
    ax.ax_row_dendrogram.set_visible(False)

    #ax = ax.ax_row_dendrogram.set_xlim([0,10])
    #ax = ax.ax_row_dendrogram.set_ylim([0,10])
    ax.ax_cbar.set_visible(False)
    ax.ax_col_dendrogram.set_visible(False)
    plt.tight_layout()

    #rcParams['figure.figsize'] = 1500, 1500

    f = ax.savefig(outfile)
    return f

def main(active_kinase_table, out_file_all, out_file_pvc, out_file_eu):

    df = pd.read_csv(active_kinase_table, sep='\t')
    domains_all = df['Signature accession merged']
    species = df['organism_name'].drop_duplicates().values.tolist()

    ## table for plotting the heatmap with species group and domain id/annotation columns
    df_domains_for_plotting = df[['group', 'organism_name', 'Signature accession merged','Signature description merged']]

    ## count the abundance of domains
    cnt_dom = count_abundance(df_domains_for_plotting)
    cnt_dom = absence_precence_table(cnt_dom)
    cnt_dom_pvc = cnt_dom[~cnt_dom['group'].isin(['Nematoda', 'Ascomycota'])]
    cnt_dom_eu = cnt_dom[cnt_dom['group'].isin(['Nematoda', 'Ascomycota'])]
    cnt_dom_matrix_df_all = cnt_dom_matrix(cnt_dom)
    cnt_dom_matrix_df_pvc = cnt_dom_matrix(cnt_dom_pvc)
    cnt_dom_matrix_df_eu = cnt_dom_matrix(cnt_dom_eu)

    f_all = create_heatmap(cnt_dom_matrix_df_all, out_file_all)
    f_pvc = create_heatmap(cnt_dom_matrix_df_pvc, out_file_pvc)
    f_eu = create_heatmap(cnt_dom_matrix_df_eu, out_file_eu)
    return f_all, f_pvc, f_eu


# get input and run the script

active_kinase_table = snakemake.input.table
out_file_all = snakemake.output.heatmap_all
out_file_pvc = snakemake.output.heatmap_pvc
out_file_eu = snakemake.output.heatmap_eu
f1, f2, f3= main(active_kinase_table, out_file_all, out_file_pvc, out_file_eu)
