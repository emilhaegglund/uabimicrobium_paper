import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import os
import pandas as pd

def open_file(my_file):
    """
    open the hmm file as a dataframe
    """

    col_names = ['target name','accession','tlen', 'query name','accession domain','qlen','E-value', 'score', 'bias','#', 'of', 'c-Evalue', 'i-Evalue', 'score2', 'bias2', 'domain from', 'domain to', 'protein from', 'protein to']

    df = pd.read_csv(my_file, sep='\s',comment='#', header=None, names= col_names, usecols=col_names, index_col=False)
    df = df.fillna('Missing')

    return df

global d_K
d_K = { "Brocadiaceae" : 30,
        "Isosphaeraceae": 33,
        "Pirellula": 54,
        "Saltatorellus": 36,
        "Bythopirellula": 46,
        "Orphan": 61,
        "Gemmataceae": 54,
        "Outgroup": 34,
        "Phycisphaerae": 29,
        "Gimesia": 31,
        "Nematoda": 87,
        "Ascomycota": 88
        }

global d_RD
d_RD = { "Brocadiaceae" : 174,
        "Isosphaeraceae": 290,
        "Pirellula": 356,
        "Saltatorellus": 255,
        "Bythopirellula": 215,
        "Orphan": 644,
        "Gemmataceae": 480,
        "Outgroup": 257,
        "Phycisphaerae": 203,
        "Gimesia": 292,
        "Nematoda": 524,
        "Ascomycota": 452
        }

global d_D1
d_D1 = { "Brocadiaceae" : 175,
        "Isosphaeraceae": 291,
        "Pirellula": 357,
        "Saltatorellus": 256,
        "Bythopirellula": 216,
        "Orphan": 654,
        "Gemmataceae": 481,
        "Outgroup": 258,
        "Phycisphaerae": 204,
        "Gimesia": 293,
        "Nematoda": 525,
        "Ascomycota": 453
        }

global d_D2
d_D2 = { "Brocadiaceae" : 200,
        "Isosphaeraceae": 372,
        "Pirellula": 445,
        "Saltatorellus": 320,
        "Bythopirellula": 275,
        "Orphan": 735,
        "Gemmataceae": 602,
        "Outgroup": 293,
        "Phycisphaerae": 234,
        "Gimesia": 343,
        "Nematoda": 605,
        "Ascomycota": 606
        }


def extract_active_sites(file_path, group):
    data = []
    for record in SeqIO.parse(file_path, "fasta"):

        protein_id = record.id.split('_')[0]
        domain_number = record.id.split('_')[-1]
        domain_id = record.id

        catalytic_site_lys = str(record.seq[d_K[group]])
        catalytic_site_asp_1 = str(record.seq[d_D1[group]])
        catalytic_site_asp_2 = str(record.seq[d_D2[group]])
        RD_site = str(record.seq[d_RD[group]])

        data.append([protein_id, domain_id, int(domain_number), catalytic_site_lys, catalytic_site_asp_1, catalytic_site_asp_2, RD_site, group])
    return data


def add_status(df_per_domain):
    status_table = []

    for row_index in df_per_domain.index.values:
        row = df_per_domain.iloc[[row_index]]

        if int(row['domain_from'].values[0]) <= 29 and int(row['domain_to'].values[0])>=140 and row['catalytic_site_lys'].values[0]=='K' and row['catalytic_site_asp_1'].values[0]=='D' and row['catalytic_site_asp_2'].values[0]=='D':
            status_table.append([row['domain_id'].values[0], 'full'])
        else:
            status_table.append([row['domain_id'].values[0], 'fraction'])

    status_df = pd.DataFrame(status_table, columns = ['domain_id', 'domain_completness'])
    df_per_domain = df_per_domain.merge(status_df, how='left', on='domain_id')

    return df_per_domain


def filter_out_inactive_fraction_goups(df_per_domain):

    df = df_per_domain[df_per_domain['domain_completness']=='fraction']
    df = df.reset_index()
    ## remove the kdd sites outside range

    for row_index in df.index:
        row = df.iloc[row_index,:]

        KDD_range = list(range(row['domain_from'], row['domain_to']))
        if 29 not in KDD_range and row['catalytic_site_lys'] == 'K':
            df.loc[row_index,'catalytic_site_lys'] = 'X'
        if 122 not in KDD_range and row['catalytic_site_asp_1'] == 'D':
            df.loc[row_index,'catalytic_site_asp_1'] = 'X'
        if 140 not in KDD_range and row['catalytic_site_asp_2'] == 'D':
            df.loc[row_index,'catalytic_site_asp_2'] = 'X'

    df = df.assign(KDD_in_range = df.domain_from.astype(str) + ', ' + df.domain_to.astype(str))

    df['catalytic_site_lys'] = df['catalytic_site_lys'].astype(str) + '1'
    df['catalytic_site_asp_1'] = df['catalytic_site_asp_1'].astype(str) + '2'
    df['catalytic_site_asp_2'] = df['catalytic_site_asp_2'].astype(str) + '3'

    df = df.assign(KDD = df.catalytic_site_lys.astype(str) + ', ' + df.catalytic_site_asp_1.astype(str) + ', ' + df.catalytic_site_asp_2.astype(str))


    df_1 = df[['protein_id', 'KDD']].groupby(['protein_id'])['KDD'].apply(','.join).reset_index()
    df_1 = df_1[df_1['KDD'].str.contains('K1')]
    df_1 = df_1[df_1['KDD'].str.contains('D2')]
    df_1 = df_1[df_1['KDD'].str.contains('D3')]

    active_frac_pid = df_1['protein_id'].values.tolist()
    df_act = df[df['protein_id'].isin(active_frac_pid)]
    active_frac_did = df_act['domain_id'].values.tolist()

    df = df[['protein_id', 'domain_id', 'domain_number', 'catalytic_site_lys',
       'catalytic_site_asp_1', 'catalytic_site_asp_2', 'RD_site', 'group',
       'domain_from', 'domain_to', 'protein_from', 'protein_to',
       'domain_completness']]
    return active_frac_pid, active_frac_did, df_act, df

def merge_hmm_files(hmm_folder, constraint):
    counter = 0
    df_hmm_all = []
    for hmm_file in os.listdir(hmm_folder):
        hmm_path = os.path.join(hmm_folder, hmm_file)
        hmm_df = open_file(hmm_path)[['target name', '#', 'domain from', 'domain to', 'protein from', 'protein to']]
        df_hmm_all.append(hmm_df)
    df_hmm_all = pd.concat(df_hmm_all)
    df_hmm_all.columns = ['protein_id', 'domain_number', 'domain_from', 'domain_to', 'protein_from', 'protein_to']
    return df_hmm_all

def extract_kdd_sites(msa_folder, groups):
    table = []

    for msa_file in os.listdir(msa_folder):
        group = msa_file.split('_')[-1].split('.aln')[0]
        if msa_file.split('_')[-1].split('.aln')[0] in groups:
            input_file = msa_folder + msa_file
            data = extract_active_sites(input_file, group)
            table+=data

    df_per_domain = pd.DataFrame(table, columns = ['protein_id', 'domain_id', 'domain_number', 'catalytic_site_lys','catalytic_site_asp_1','catalytic_site_asp_2','RD_site', 'group'])

    return df_per_domain

def main(msa_folder, groups, hmm_folder, output_1, output_2):

    # open the hmm and msa files and merge them for all species/groups
    df_hmm_all = merge_hmm_files(hmm_folder, 'reformated.tsv')
    df_per_domain = extract_kdd_sites(msa_folder, groups)

    # merge the hmm df and the KDD info df
    df_per_domain = df_per_domain.merge(df_hmm_all, how='left', on=['protein_id','domain_number'])

    # add the status of the domain (complete or fraciton)
    df_per_domain = add_status(df_per_domain)

    # filter out the non active protein kinases and save the result as a table
    pid_active_frac, did_active_frac, df_act, df = filter_out_inactive_fraction_goups(df_per_domain)
    did_active = did_active_frac + df_per_domain[df_per_domain['domain_completness']=='full']['domain_id'].drop_duplicates().values.tolist()

    #merge active fractions and remove fractions that does not build up to a whole domain.
    #df_per_domain = df_per_domain[df_per_domain['domain_id'].isin(did_active)]


    # save a table with active and non active
    df_per_domain_active = df_per_domain[df_per_domain['domain_id'].isin(df_per_domain[df_per_domain['domain_completness']=='full']['domain_id'].drop_duplicates().values.tolist())]
    df_per_domain_active['Function'] = 'active'
    df_per_domain_possibly_active = df_per_domain[df_per_domain['domain_id'].isin(did_active_frac)]
    df_per_domain_possibly_active['Function'] = 'possibly active'

    df_per_domain_non_active = df[~df['domain_id'].isin(did_active)]
    df_per_domain_non_active['Function'] = 'non-active'
    df_per_domain_all = pd.concat([df_per_domain_active,df_per_domain_possibly_active, df_per_domain_non_active],ignore_index=True)

    df_per_domain.to_csv(output_1, index=False, sep='\t')
    df_per_domain_all.to_csv(output_2, index=False, sep='\t')

    return df_per_domain

msa_folder = snakemake.params.msa_folder
group_list = snakemake.params.group_list
hmm_folder = snakemake.input.hmm_folder
output_1 = snakemake.output.output_1
output_2 = snakemake.output.output_2

df_domain = main(msa_folder, group_list, hmm_folder, output_1, output_2)
