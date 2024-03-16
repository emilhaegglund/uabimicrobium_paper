import pandas as pd

def find_domains_to_merge(df):
    pids = df['protein_id'].drop_duplicates().values.tolist()
    merge_list = []
    pid_merged = []

    for pid in pids:
        pid_df = df[df['protein_id']==pid].reset_index()
        row_index_list = pid_df.index.values.astype(int)
        first = 0

        for row_index in row_index_list:
            domain_df = pid_df.iloc[row_index]

            if first == 0:
                before_id = domain_df['domain_id']
                before_ds = domain_df['domain_from']
                before_de = domain_df['domain_to']
                before_ps = domain_df['protein_from']
                before_pe = domain_df['protein_to']
                before_completness = domain_df['domain_completness']
                first = 1

            else:
                current_id = domain_df['domain_id']
                current_ds = domain_df['domain_from']
                current_de = domain_df['domain_to']
                current_ps = domain_df['protein_from']
                current_pe = domain_df['protein_to']
                current_completness = domain_df['domain_completness']

                if int(current_de) > int(before_de) and int(current_ds) > int(before_ds):

                    #if current_completness != before_completness or (current_completness == 'fraction' and before_completness=='fraction'):
                    if current_completness == before_completness and current_completness == 'fraction':
                        merge_list.append([pid, before_id,current_id, before_ds, current_de,before_ps, current_pe])

                        if pid not in pid_merged:
                            pid_merged.append(pid)

                before_id = domain_df['domain_id']
                before_ds = domain_df['domain_from']
                before_de = domain_df['domain_to']
                before_ps = domain_df['protein_from']
                before_pe = domain_df['protein_to']
                before_completness = domain_df['domain_completness']

    merge_domain_df = pd.DataFrame(merge_list, columns = ['protein_id','domain_1','domain_2', 'domain_from', 'domain_to', 'protein_from', 'protein_to'])
    domain_id_fraction_to_keep = merge_domain_df['domain_1'].drop_duplicates().values.tolist() + merge_domain_df['domain_2'].drop_duplicates().values.tolist()
    domain_id_full_to_keep = df[df['domain_completness']=='full']['domain_id']
    domain_id_to_keep = domain_id_full_to_keep.drop_duplicates().values.tolist() + domain_id_fraction_to_keep
    df_new = df[df['domain_id'].isin(domain_id_to_keep)]

    return merge_domain_df, df_new


def merge_domains(merge_domain_df):
    """
    merge the fractional domains and give them new domain ids and start, stop values.
    """
    pids = merge_domain_df['protein_id'].drop_duplicates().values.tolist()
    merge = []

    for pid in pids:
        first_dom='yes'
        pid_df = merge_domain_df[merge_domain_df['protein_id']==pid].reset_index()
        pid_index_list = pid_df.index.values.astype(int)

        for row_index in pid_index_list:

            dom_to_merge = pid_df.iloc[row_index]
            domain_id_new = str(dom_to_merge['domain_1']) + ';' + str(dom_to_merge['domain_2'])
            d_start = dom_to_merge['domain_from']
            d_end = dom_to_merge['domain_to']
            p_start = dom_to_merge['protein_from']
            p_end = dom_to_merge['protein_to']

            # if the first domain merge the two specified domains
            if first_dom == 'yes':
                merge.append([pid, domain_id_new, d_start, d_end, p_start, p_end, dom_to_merge['domain_2']])
                first_dom = 'No'

            else:

                # check if a third or more fraction will be merged
                if dom_to_merge['domain_1'] == merge[-1][-1]:
                    before = merge[-1]
                    merge.append([pid, str(before[1])+ ';' + str(dom_to_merge['domain_2']), before[2], d_end, before[4], p_end, dom_to_merge['domain_2']])
                # if only two domains will be merged
                else:
                    merge.append([pid, domain_id_new, d_start, d_end, p_start, p_end, dom_to_merge['domain_2']])

    df_merged_fractional_domains = pd.DataFrame(merge, columns = ['protein_id', 'domain_id_new', 'domain_from', 'domain_to', 'protein_from', 'protein_to', 'last domain'])
    return df_merged_fractional_domains


def KDD_column(df):

    df['catalytic_site_lys'] = df['catalytic_site_lys'].astype(str) + '1'
    df['catalytic_site_asp_1'] = df['catalytic_site_asp_1'].astype(str) + '2'
    df['catalytic_site_asp_2'] = df['catalytic_site_asp_2'].astype(str) + '3'

    df = df.assign(KDD = df.catalytic_site_lys.astype(str) + ', ' + df.catalytic_site_asp_1.astype(str) + ', ' + df.catalytic_site_asp_2.astype(str))

    return df

def extract_kdd_from_merged(df, df_merged_fractional_domains):
    """
    return the merged domains ids of merged domains with the KDD sites
    """
    rows = df_merged_fractional_domains.index
    merged_domains_w_kdd = []

    for row in rows:
        domain_ids_list = df_merged_fractional_domains.iloc[[row]]['domain_id_new'].values.tolist()[0]
        domain_ids = domain_ids_list.split(';')
        kdd_list = []

        for domain_id in domain_ids:
            KDD = df[df['domain_id']==domain_id][['domain_id','catalytic_site_lys', 'catalytic_site_asp_1', 'catalytic_site_asp_2']]
            KDD = KDD_column(KDD)['KDD'].values.tolist()
            kdd_list += KDD[0].split(', ')

        if 'K1' in kdd_list and 'D2' in kdd_list and 'D3' in kdd_list:
            merged_domains_w_kdd.append(domain_ids_list)

    return merged_domains_w_kdd

def get_rd_protein_ids(df):
    pids = df[df['RD_site']=='R']['protein_id'].drop_duplicates().values.tolist()
    return pids

def append_info(df, df_merged):
    #df_merged = df_merged_fractional_domains[df_merged_fractional_domains['domain_id_new'].isin(merged_domains_w_kdd)].reset_index(drop=True)
    df_groups = df[['protein_id','group']].drop_duplicates()
    df_merged = df_merged[['protein_id','domain_id_new', 'domain_from', 'domain_to', 'protein_from', 'protein_to']]
    df_merged.columns = ['protein_id','domain_id', 'domain_from', 'domain_to', 'protein_from', 'protein_to']
    df_merged['catalytic_site_lys'] = 'K'
    df_merged['catalytic_site_asp_1'] = 'D'
    df_merged['catalytic_site_asp_2'] = 'D'

    RD_site = []
    for pid in df_merged['protein_id'].drop_duplicates().values.tolist():
        site = df[df['protein_id']==pid]['RD_site'].drop_duplicates().values.tolist()
        site = ';'.join(site)
        RD_site.append([pid, site])
    RD_site_df = pd.DataFrame(RD_site)

    RD_site_df.columns = ['protein_id', 'RD_site']

    df_merged = df_merged.merge(RD_site_df, how='left',  on='protein_id')
    df_merged['domain_completness'] = 'full merged'
    df_merged = df_merged.merge(df_groups, how='left',  on='protein_id')
    return df_merged

def main(table_path, outfile_name):
    df = pd.read_csv(table_path, sep='\t')
    merge_domain_df, df_to_keep = find_domains_to_merge(df)
    df_merged_fractional_domains = merge_domains(merge_domain_df)

    # remove merged fractional domains based on length (missing K or last d site)
    df_merged_fractional_domains = df_merged_fractional_domains[df_merged_fractional_domains['domain_from']<=29]
    df_merged_fractional_domains = df_merged_fractional_domains[df_merged_fractional_domains['domain_to']>=140]

    # remove merged fractional domains if they do not have the KDD sites
    df_merged_fractional_domains = df_merged_fractional_domains.reset_index(drop=True)
    merged_domains_w_kdd = extract_kdd_from_merged(df, df_merged_fractional_domains)
    df_merged = df_merged_fractional_domains[df_merged_fractional_domains['domain_id_new'].isin(merged_domains_w_kdd)].reset_index(drop=True)

    # add info to the merged domains and merge with the full domains table
    df_merged = append_info(df, df_merged)
    df_full = df[df['domain_completness']=='full'][['protein_id','domain_id','catalytic_site_lys','catalytic_site_asp_1','catalytic_site_asp_2','RD_site', 'group','domain_from','domain_to','protein_from','protein_to', 'domain_completness']]
    #df_full = df_full[df_full['domain_id'].isin(merged_domains_w_kdd)]
    df_new = pd.concat([df_merged,df_full])

    # save to file
    f = df_new.to_csv(outfile_name, sep='\t', index=False)
    return df_new

df_new = main(snakemake.input[0], snakemake.output[0])


