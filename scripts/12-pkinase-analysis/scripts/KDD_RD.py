import pandas as pd
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
import seaborn as sns

def open_files(phobius, active_sites_table, p_mapp, genome_stat, genome_stat_eu):
    df_tm = pd.read_csv(phobius, sep='\t')
    print(df_tm.columns)
    df_info = df_tm[['species','nr of pkinases', 'all proteins']]
    df_tm = df_tm[['species', 'group','% transmembrane proteins', 'transmembrane proteins']]
    df_info['% kinases'] = round((df_info['nr of pkinases']/df_info['all proteins'])*100, 1)
    df_rd = pd.read_csv(active_sites_table, sep='\t')
    print(df_rd.columns)
    df_pmapp = pd.read_csv(p_mapp, sep='\t')
    df_genome =  pd.read_csv(genome_stat, sep='\t')[['assembly_accession','organism_name']]
    df_eu_genomes = pd.read_csv(genome_stat_eu, sep='\t')[['accession','organism name', 'protein_id']]
    df_eu_genomes.columns = ['accession', 'species', 'protein_id']

    return df_tm, df_rd, df_pmapp, df_genome, df_info, df_eu_genomes


def add_species(df_rd, df_pmapp, df_genome, df_eu_genomes):
    
    df_genome.columns = ['accession', 'species']
    df_genome = df_genome.merge(df_pmapp, how='left', on = 'accession')
    df_genome = pd.concat([df_genome, df_eu_genomes], ignore_index=True)
    df_rd = df_rd.merge(df_genome, how='left', on = 'protein_id')

    return df_rd

def change_uab_group_name(df_rd):
    uab_pid = df_rd[df_rd['protein_id'].str.contains('BBM')]['protein_id'].drop_duplicates().values.tolist()
    f = lambda x: 'Orphan_Uab' if (x in uab_pid) else df_rd[df_rd['protein_id']==x]['group'].values[0]
    df_rd['group'] = df_rd['protein_id'].map(f)
    return df_rd

def count_rd_proteins(df_rd,df_info):
    df_info = df_info[['species', 'nr of pkinases']]
    df_rd = df_rd[['protein_id', 'RD_site', 'species', 'group']]
    df_rd = df_rd[df_rd['RD_site'].str.contains('R')]
    df_rd = df_rd[['protein_id','species', 'group']].drop_duplicates()
    df_rd = df_rd[['species','group']].value_counts().to_frame().reset_index()
    df_rd.columns = ['species','group','RD proteins']
    
    df_rd = df_rd.merge(df_info, how='left', on='species')
    df_rd['% RD proteins'] = (df_rd['RD proteins'] / df_rd['nr of pkinases'])*100
    df_rd['% RD proteins'] = df_rd['% RD proteins'].round(0)
    df_rd = df_rd[['species','group','% RD proteins', 'RD proteins']]
    
    return df_rd

def merge_rd_tm(df_rd, df_tm, df_info):
    df_rd['type'] = 'RD'
    df_rd.columns = ['species','group', 'protein count frac', 'protein count','type']
    df_tm['type'] = 'TM'
    df_tm.columns = ['species','group', 'protein count frac', 'protein count', 'type']
    df_tm_rd = pd.concat([df_rd, df_tm], ignore_index=True)
    
    df_no_rd = df_tm[df_tm['protein count']==0].drop_duplicates()
    df_no_rd['type'] = 'RD'
    df_tm_rd = pd.concat([df_tm_rd,df_no_rd], ignore_index=True)
    
    df = df_tm_rd.merge(df_info, how='left', on='species')
    df['all proteins'] = df['all proteins'].round(0)
    df['% kinases'] = df['% kinases'].round(0)
    
    df = df.astype({'all proteins':'int', '% kinases':'int'})
    return df

def barplot(df, tm_figure_name):
    
    group_order = ['Orphan_Uab', 'Saltatorellus', 'Brocadiaceae', 'Phycisphaerae','Gemmataceae', 'Isosphaeraceae', 'Gimesia', 'Bythopirellula', 'Pirellula', 'Orphan', 'Outgroup', 'Nematoda', 'Ascomycota', 'Proteobacteria']
    #df = df.sort_values(by="group", key=lambda column: column.map(lambda e: group_order.index(e)))
    
    df_test = df
    df_test = df_test[df_test['type']=='TM']
    df_group_tm_mean = df_test.groupby(by='group')['protein count'].mean()
    df_group_tm_mean = df_group_tm_mean.reset_index()
    df_group_tm_mean.columns = ['group','protein count mean']
    
    df_test = df_test.merge(df_group_tm_mean, how='left', on='group')
    df_sorted_groups = df_test.sort_values(by='protein count mean', ascending=False)
    df_sorted_groups = df_sorted_groups['group'].drop_duplicates().values.tolist()

    df = df.sort_values(by="group", key=lambda column: column.map(lambda e: df_sorted_groups.index(e)))
    
    plt.rcParams["figure.figsize"] = [8.3,5.3]
    plt.rcParams["figure.autolayout"] = True

    ax = sns.barplot(data=df, x="protein count", y="group", hue="type", palette=['#C4D7EB', '#6B9CCE'])
    #sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    #ax.set_yscale('log')
    #plt.xticks(rotation=90)
    f = ax.get_figure() 
    #figure(figsize=(8.3,11.7), dpi=80)
    #f.ax.legend(loc=2)
    plt.legend(bbox_to_anchor=(1.25, 1), borderaxespad=0)
    f.savefig(tm_figure_name)
    return f
    
#def add_non_tm_species(df, df_genome, df_eu_genomes):
#    species_present = df.species.drop_duplicates().values.tolist()
#    species_abscent_pvc = df_genome[~df_genome['assembly_accession'].isin(species_present)]
#    species_abscent_eu = df_eu_genomes[~df_eu_genomes['accession'].isin(species_present)]
    
#    return df


def create_df_info(protein_mapping, genome_stat, genome_stat_p_mapp_eu):
    p_mapp_df = pd.read_csv(protein_mapping, sep='\t')
    genome_stat_df = pd.read_csv(genome_stat, sep='\t')[['organism_name', 'group', 'assembly_accession']]
    genome_stat_df.columns = ['organism_name', 'group', 'accession']
    genome_stat_p_mapp_eu_df = pd.read_csv(genome_stat_p_mapp_eu, sep='\t')[['organism name', 'group', 'accession', 'protein_id']]
    genome_stat_p_mapp_eu_df.columns = ['organism_name', 'group', 'accession', 'protein_id']
    
    df_info = genome_stat_df.merge(p_mapp_df, how='left', on='accession')
    df_info = pd.concat([df_info,genome_stat_p_mapp_eu_df])[['organism_name', 'accession', 'group', 'protein_id']]
    return df_info

def main(phobius, active_sites_table, p_mapp, genome_stat, genome_stat_eu, kdd_rd_table, tm_figure_name):
    df_tm, df_rd, df_pmapp, df_genome, df_info,df_eu_genomes = open_files(phobius, active_sites_table, p_mapp, genome_stat, genome_stat_eu)
    df_rd = add_species(df_rd, df_pmapp, df_genome, df_eu_genomes)
    
    df_rd = change_uab_group_name(df_rd)
    
    df_rd = count_rd_proteins(df_rd, df_info)
    df  = merge_rd_tm(df_rd, df_tm, df_info)
    
    #df = add_non_tm_species(df, df_genome, df_eu_genomes)
    
    df.to_csv(kdd_rd_table, sep='\t', index=False)
    f = barplot(df, tm_figure_name)
    
    
    return df, f

df = main(snakemake.input.phobius, snakemake.input.active_sites_table, snakemake.input.p_mapp, snakemake.input.genome_stat, snakemake.input.genome_stat_eu, snakemake.output.kdd_rd_table, snakemake.output.tm_figure_name)
