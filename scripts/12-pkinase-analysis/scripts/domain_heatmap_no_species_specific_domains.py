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
    cnt_dom.columns = ['group','organism_name', 'Signature accession merged', 'Signature description merged', 'Count']
    cnt_dom = cnt_dom[['group','organism_name','Signature description merged', 'Count']]
    #cnt_dom = cnt_dom[cnt_dom['Count']>1.0]
    return cnt_dom

def absence_precence_table(cnt_dom):
    cnt_abs_pres = cnt_dom
    cnt_abs_pres['Count'].values[cnt_abs_pres['Count'] > 1.0] = 1.0
    cnt_abs_pres['Count'].values[cnt_abs_pres['Signature description merged'] == 'Missing'] = 0
    cnt_abs_pres.drop_duplicates()

    return cnt_abs_pres

def cnt_dom_matrix(cnt_df):
    
    sort_df = cnt_df
    cnt_dom_matrix = cnt_df.pivot_table(index=['organism_name', 'group'], columns='Signature description merged',values='Count',sort=False).fillna(0)
    cnt_dom_matrix = cnt_dom_matrix.T
    #print(cnt_dom_matrix)
    cnt_dom_matrix = sort_heatmap(cnt_dom_matrix, sort_df)
    
    return cnt_dom_matrix

def sort_heatmap(df, cnt_df):
    
    cnt_df_1 = cnt_df.copy(deep=True)
    cnt_df_2 = cnt_df.copy(deep=True)
    
    # sort the species order
    species_order = ['Planctomycetes bacterium SRT547','Planctomycetes bacterium Poly30', 'Planctomycetes bacterium Pla163', 'Planctomycetes bacterium Pla133', 'Planctomycetes bacterium Pla86', 'Candidatus Scalindua japonica', 'Candidatus Brocadiaceae bacterium S225', 'Candidatus Kuenenia stuttgartiensis','Candidatus Jettenia caeni','Candidatus Brocadia sinica JPN1','Candidatus Brocadia sapporoensis', 'Candidatus Brocadiaceae bacterium B188', 'Anaerohalosphaera lusitana', 'Sedimentisphaera salicampi', 'Phycisphaerae bacterium RAS2',  'Phycisphaera mikurensis NBRC 102666', 'Planctomycetes bacterium Pan265', 'Tuwongella immobilis', 'Limnoglobus roseus', 'Urbifossiella limnaea', 'Gemmata sp. SH-PL17', 'Gemmata obscuriglobus', 'Planctomycetes bacterium Pan216','Isosphaera pallida ATCC 43644', 'Tautonia plasticadhaerens','Singulisphaera acidiphila DSM 18658', 'Planctomyces sp. SH-PL62','Aquisphaera giovannonii', 'Symmachiella dynata', 'Planctomycetes bacterium Pan189', 'Caulifigura coniformis','Gimesia chilikensis', 'Gimesia maris', 'Planctomycetales bacterium 10988', 'Thermogutta terrifontis', 'Bythopirellula goksoyri', 'Lacipirellula parvula', 'Pirellulimonas nuda', 'Planctomycetes bacterium MalM25', 'Botrimarina mediterranea', 'Bremerella volcania', 'Mariniblastus fucicola', 'Rhodopirellula baltica SH 1', 'Stieleria maiorica','Planctomycetes bacterium FF011L','Victivallales bacterium CCUG 44730', 'Akkermansia muciniphila', 'Opitutus terrae PB90-1', 'Nibricoccus aquaticus', 'Akkermansia muciniphila ATCC BAA-835', 'Coraliomargarita akajimensis DSM 45221', 'Candidatus Protochlamydia naegleriophila', 'Simkania negevensis Z', 'Waddlia chondrophila WSU 86-1044', 'Akkermansia glycaniphila', 'Neochlamydia sp. S13', 'Lacunisphaera limnophila', 'Opitutaceae bacterium TAV5', 'Parachlamydia acanthamoebae UV-7', 'Kiritimatiella glycovorans', 'Verrucomicrobia bacterium S94', 'Verrucomicrobia bacterium IMCC26134', 'Methylacidiphilum infernorum V4','Methylacidiphilum kamchatkense Kam1','Candidatus Xiphinematobacter sp. Idaho Grape', 'Caenorhabditis elegans', 'Saccharomyces cerevisiae', 'Schizosaccharomyces pombe']
    
    df = df.sort_index(axis=1, level='organism_name', key=lambda column: column.map(lambda e: species_order.index(e)))
    
    # sort domain order all
    domain_order_sp = cnt_df[['Signature description merged', 'Count']]
    domain_order_sp = domain_order_sp.groupby('Signature description merged')['Count'].sum()
    domain_order_sp = domain_order_sp.to_frame().reset_index()
    domain_order_sp.columns = ['Signature description merged','Count sp']
    ##domain_order = domain_order.sort_values(by='Count', ascending=False)
    domain_order = cnt_df[['Signature description merged']].value_counts().reset_index()
    domain_order.columns = ['Signature description merged','Count']
    domain_order = domain_order.merge(domain_order_sp, how='left', on='Signature description merged')
    domain_order = domain_order.sort_values(by=['Count sp','Count'], ascending=False)
    domain_order = domain_order['Signature description merged'].drop_duplicates().values.tolist()
    df = df.sort_values(by="Signature description merged", key=lambda column: column.map(lambda e: domain_order.index(e)))

    # sort the domain order species/group specific
    df_temp = df.copy(deep=True)
    single_domains = cnt_df_2[['Signature description merged', 'Count']]['Signature description merged'].drop_duplicates(keep=False).values.tolist()
    
    single_domain_order = []
    single_domain_order_temp = []
    species = df.columns.get_level_values(0).values.tolist()
    col_index = 0
    
    for specie in species:
        col = df_temp.iloc[:,col_index]
        single_domain_order_temp=[]
        if col.empty == False:
            ind = col.index.values.tolist()
            for i in ind:
                if i in single_domains:
                    dom = col[col.index.str.contains(i, regex=False)]
                    if dom.empty == False:
                        if dom.values.tolist()[0] > 0:
                            single_domain_order_temp.append([i,dom.values.tolist()[0]])
        sdo = pd.DataFrame(single_domain_order_temp, columns=['dom_index','count']).sort_values(by='count', ascending=False)['dom_index'].values.tolist()
        single_domain_order+=sdo
        col_index+=1
        
    domains = df.index.values.tolist() 
    domains_other = []

    for dom in domains:
        if dom not in single_domain_order:
            domains_other.append(dom)

    single_domain_order = domains_other + single_domain_order
    
    df = df.sort_values(by="Signature description merged", key=lambda column: column.map(lambda e: single_domain_order.index(e)))
    


    return df

def create_heatmap(cnt_dom_matrix, outfile):


    cnt_dom_matrix = cnt_dom_matrix.replace(0, np.nan)
    mask = cnt_dom_matrix.isnull()

    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    
    #c = matplotlib.colors.LinearSegmentedColormap.from_list("", ['white', '#6B9CCE']) # '#C4D7EB', '#6B9CCE'
    #c_2 = matplotlib.colors.LinearSegmentedColormap.from_list("", [ '#FFBF65','#FF96C5','#C05780','#E77577', '#6C88C4', '#00A5E3'])##FF828B']) '#00B0BA', '#4DD091', v3 = [ '#FFBF65','#FF96C5','#C05780','#E77577', '#6C88C4', '#00A5E3']
    #c_3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#d24e01', '#ca5cdd'])# '#C696C6', '#301934']) '#6B9CCE', '
    #c_4 = sns.cubehelix_palette(start=2,rot=.6,as_cmap=True) #start=-.6,rot=.6
    #c_5 = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#C4D7EB',"#6B9CCE",'#bc96ff', '#FF96C5', '#C05780'])
    #c_6 = matplotlib.colors.LinearSegmentedColormap.from_list("", ['lightyellow', 'gold', 'lightskyblue', 'royalblue'])
    c_7 = matplotlib.colors.LinearSegmentedColormap.from_list("", ['lightyellow', 'gold', 'royalblue'])
    ax = sns.heatmap(cnt_dom_matrix, cmap='viridis_r', cbar=True, linewidth=.5, linecolor='lightgray', mask=mask, annot=False) #cmap =c
    figure = ax.get_figure()
    plt.tight_layout()
    
    plt.figure(figsize=(20,35), dpi=80) # 8.3,11.7

    f = figure.savefig(outfile)
    return f

def main(genomes, active_kinase_table, out_file_all, out_file_pvc, out_file_eu):
    
    df = pd.read_csv(active_kinase_table, sep='\t')
    domains_all = df['Signature accession merged']
    species = df['organism_name'].drop_duplicates().values.tolist()

    ## table for plotting the heatmap with species group and domain id/annotation columns
    df_domains_for_plotting = df[['group', 'organism_name', 'Signature accession merged','Signature description merged']]
    
    ## Add species without protein kinases
    df_genomes = pd.read_csv(genomes, sep='\t')
    df_no_kinase_species = df_genomes[~df_genomes['organism_name'].isin(species)][['group', 'organism_name']]
    df_no_kinase_species['Signature accession merged'] = 'Missing'#np.NaN
    df_no_kinase_species['Signature description merged'] = 'Missing'#np.NaN
    print(df_no_kinase_species)
    df_domains_for_plotting = pd.concat([df_domains_for_plotting, df_no_kinase_species], ignore_index=True)
    print(df_domains_for_plotting)
    ## count the abundance of domains
    cnt_dom = count_abundance(df_domains_for_plotting)
    #cnt_dom = absence_precence_table(cnt_dom)
    
    cnt_dom_pvc = cnt_dom[~cnt_dom['group'].isin(['Nematoda', 'Ascomycota', 'Outgroup'])]
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
genomes = snakemake.input.genomes
out_file_all = snakemake.output.heatmap_all
out_file_pvc = snakemake.output.heatmap_pvc
out_file_eu = snakemake.output.heatmap_eu
f1, f2, f3= main(genomes, active_kinase_table, out_file_all, out_file_pvc, out_file_eu)