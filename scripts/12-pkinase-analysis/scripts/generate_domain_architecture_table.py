import argparse
import pandas as pd
parser = argparse.ArgumentParser()
from collections import Counter


def open_file(my_file, header, seperator, col_names):
    if col_names is not None:
        df = pd.read_csv(my_file, sep=seperator, header=header, names=col_names)
    else:
        df = pd.read_csv(my_file, sep=seperator, header=header)
    df = df.fillna('Missing')

    return df

def open_interproscan_file(interproscan_file, analysis, interproscan_header):
    """
    extract the Pfam domains of a interproscan file
    """
    if interproscan_header == []:
        interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession',
                           'Signature description','Start location', 'Stop location', 'Score', 'Status',
                           'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
                           'optional: GO terms', 'optional: Pathway annotation']

    domain_df = pd.read_csv(interproscan_file, sep='\t', names = interproscan_header)
    domain_df = domain_df[domain_df['Analysis']==analysis]
    domain_df= domain_df.astype({'Stop location':int, 'Start location': int})
    domain_df = domain_df.sort_values(by=['protein_id','Start location']).reset_index()
    domain_df = domain_df[['protein_id','Start location', 'Stop location','Signature accession', 'Signature description', 'Sequence length']]

    return domain_df

def reshape_interproscan_file(interpro_df):
    """
    Group the domains and the annotation based on the protein ids, like this:
    ProteinID     domain_id_1;domain_id_2     domain_annotation_1;domain_annotation_2
    """

    # turn the columns into text
    for column_name in ['Start location','Stop location','Signature description','Signature accession']:
        interpro_df[column_name] = interpro_df[column_name].astype(str)

    interpro_df_pid = interpro_df.groupby('protein_id').agg({'Signature accession': lambda x: ';'.join(x)})
    interpro_df_pid = interpro_df_pid.reset_index()

    interpro_df_anno = interpro_df.groupby('protein_id').agg({'Signature description': lambda x: ';'.join(x)})
    interpro_df_anno = interpro_df_anno.reset_index()

    interpro_df = interpro_df_pid.join(interpro_df_anno.set_index('protein_id'), on='protein_id')

    return interpro_df

def get_domain_arc(interpro_df):
    """
    Extract the domain arkitectures and count how often we have a domain combination
    """
    domains_all = interpro_df['Signature accession']
    cnt_domains = Counter(domains_all.values.tolist())
    domains = domains_all.drop_duplicates().values.tolist()
    return domains, cnt_domains

def create_schematic_table(domain_arc_list, cnt_domain_arc, interpro_df, interpro_original_df): #, cnt_homolog):

    schematic_df = []
    prot_id = len(domain_arc_list)
    domain_anno = []

    for domain_arc in domain_arc_list:
        amount_arc = cnt_domain_arc[domain_arc]
        protein_ids = interpro_df[interpro_df['Signature accession']==domain_arc]['protein_id'].sort_values().str.cat(sep=';').rstrip()

        if ';' in domain_arc:

            prot_domains = domain_arc.split(';')
            shorter_before = 'no'
            smaller_region = 0
            start = 0
            stop = 0

            for p in prot_domains:
                if p in ['SIGNAL', 'n_TM_c', 'c_TM_n']:
                    smaller_region +=1

            prot_len = ((len(prot_domains)-(smaller_region))*100) + (smaller_region*50) + (25*(len(prot_domains)+1))
            first = 'yes'
            for domain in prot_domains:


                # get domain annotation and start and stop values for the domain or Transmembrane/Signal regions
                if 'SIGNAL' not in domain and 'TM' not in domain:
                    if 'X' in domain:
                        domain_anno = 'potential_domain'
                    else:
                        domain_anno = interpro_original_df[interpro_original_df['Signature accession']==domain]['Signature description'].drop_duplicates().values.tolist()[0]

                    if first == 'yes':
                        start = 25
                        stop = 125
                    if first== 'no':
                        if shorter_before == 'yes':
                            start += 75
                            stop += 125
                        else:
                            start += 125
                            stop += 125
                    shorter_before = 'no'
                else:
                    if 'SIGNAL' in domain:
                        domain_anno = 'signal'
                    if 'TM' in domain:
                        domain_anno = 'transmembrane region'

                    if first=='yes':
                        start = 25
                        stop = 75
                    if first=='no':
                        if shorter_before == 'no':
                            start += 125
                            stop += 75
                        else:
                            start += 75
                            stop += 75
                    shorter_before = 'yes'


                schematic_df.append(['Architecture ' + str(prot_id), prot_len, domain, start, stop, amount_arc, protein_ids, domain_anno])
                first = 'no'
        else:
            prot_len = 150
            start = 25
            stop = 125

            #if 'SIGNAL' in domain_arc:
            #    domain_anno == 'signal'
            #if 'TM' in domain_arc:
            #    domain_anno = 'transmembrane region'
            #if 'X' in domain:
            #    domain_anno = 'potential_domain'
            #if 'SIGNAL' not in domain_arc and 'TM' not in domain_arc and 'X' not in domain_arc:
            domain_anno = interpro_original_df[interpro_original_df['Signature accession']==domain_arc]['Signature description'].drop_duplicates().values.tolist()[0]
            schematic_df.append(['Architecture ' + str(prot_id), prot_len, domain_arc, start, stop, amount_arc, protein_ids, domain_anno])

        prot_id-=1 # the protein id counter must descrese here to get the correct order in the figure, smallest number at the top of the figure
    schematic_df = pd.DataFrame(schematic_df)

    return schematic_df

def number_of_dom(domain_ark):
    """
    count the number of domains in the domain arkitecture (domain_id_1;domain_id:2)
    """
    val = 1
    if ';' in domain_ark:
        val = len(domain_ark.split(';'))
    return val

def dom_loc(domain):
    """
    return the postion of the domain we want to sort the table on.
    """
    val = 1
    if ';' in domain:
        domain = domain.split(';')
        for e in domain:
            if e in selected_domains:
                val +=1
                return val
            val+=1
    return val


def add_potential_domains(interpro_df, aa_threshold):
    """
    Add a diamond shape if there is more than x amino acids (aa_threshold) between domain regions
    """
    potential_domain = []

    for protein in interpro_df['protein_id'].drop_duplicates():

        p = interpro_df[interpro_df['protein_id']==protein]
        start = p['Start location'].values.tolist()
        stop = p['Stop location'].values.tolist()

        p_length =interpro_df[interpro_df['protein_id']==p['protein_id'].values[0]]['Sequence length']
        i = 0
        region_length = 0
        first = 'y'

        for i in range(len(p)):

            if first == 'y':

                region_length = int(start[i]) - 0
                if region_length >= aa_threshold:
                    potential_domain.append([p['protein_id'].values[0], 0, start[i], 'X', 'Potential domain'])

            if first == 'n':
                if i == len(p):
                    region_length = int(p_length.values[0]) - int(stop[i-1])
                if i != len(p):
                    region_length = int(start[i]) - int(stop[i-1])

                if region_length >= aa_threshold:
                    potential_domain.append([p['protein_id'].values[0],stop[i-1],start[i], 'X', 'Potential domain'])

            first = 'n'

    potential_domain_df = pd.DataFrame(potential_domain, columns = ['protein_id','Start location', 'Stop location', 'Signature accession', 'Signature description'])
    interpro_df = pd.concat([interpro_df, potential_domain_df], ignore_index=True)
    interpro_df = interpro_df.sort_values(by=['protein_id','Start location']).reset_index()
    return interpro_df

def merge_similar_domains(domains_to_merge, interpro_df):
    """
    replace the chosen domain ids in the interposcan file with those in domain_to_merge file.
    """

    # open the file with domains to merge
    domains_to_merge_df = pd.read_csv(domains_to_merge, sep='\t')
    domain_id_group_df = domains_to_merge_df[['Signature accession','Signature description', 'Signature accession merged', 'Signature description merged']]

    # get the domain ids that will be changed from the domains_to_merge file
    acc_list = domain_id_group_df['Signature accession'].values.tolist()
    desc_list = domain_id_group_df['Signature description'].values.tolist()

    # extract the merged domain id and merged domain annotation that we want to use
    acc_new_list = domain_id_group_df['Signature accession merged'].values.tolist()
    anno_new_list = domain_id_group_df['Signature description merged'].values.tolist()

    # replace the chosen domain ids with the merged domain ids, chosen domain ids in acc_list
    interpro_df["Signature accession"] = interpro_df["Signature accession"].replace(acc_list,acc_new_list)
    interpro_df["Signature description"] = interpro_df["Signature description"].replace(desc_list,anno_new_list)

    return interpro_df

def add_tm_signal(phobius_file, interpro_df, pid_list):
    """
    """

    phobius_df = pd.read_csv(phobius_file, sep='\t')
    phobius_df.columns = ['protein_id','Start location', 'Stop location', 'Signature description'] ##### Switch
    phobius_df = phobius_df[phobius_df['protein_id'].isin(pid_list)]

    # change here to correct description?
    phobius_df['Signature accession'] = phobius_df['Signature description']

    phobius_df['Sequence length'] = phobius_df['Stop location'] - phobius_df['Start location']

    interpro_tm_df = pd.concat([interpro_df, phobius_df], ignore_index=True)
    interpro_tm_df = interpro_tm_df.sort_values(by=['protein_id','Start location']).reset_index()
    return interpro_tm_df


## MAIN

def main(interproscan_file, analysis, out_file, merged_domains, phobius, pid_list, groups, interproscan_header, schematic, potential_domains, sort_on_domain):

    # open interproscan df (if you do not specify the interproscan header a default header is used)
    interpro_original_df = open_interproscan_file(interproscan_file, analysis, interproscan_header)
    interpro_df = interpro_original_df

    ### Which proteins you wish to include from the interproscan_file (OPTIONAL)
    if pid_list != []:
        interpro_df = interpro_df[interpro_df['protein_id'].isin(pid_list)]


    ### merge domain ids that are the same domain (OPTIONAL)
    if merged_domains != []:
        interpro_df = merge_similar_domains(merged_domains, interpro_df)
        interpro_original_df = interpro_df

    ### add tm and signaling regions (OPTIONAL)
    if phobius != []:
        interpro_df = add_tm_signal(phobius, interpro_df, pid_list)

    ### add potential domains (OPTIONAL if schematic domains is choosen)
    #if schematic == 'yes':
    #    if potential_domains != 0:
    #        interpro_df = add_potential_domains(interpro_df, potential_domains)

    # check overlapping domains / regions
    over_lapp = interpro_df[['protein_id', 'Start location', 'Stop location']]

    for protein_id in over_lapp.protein_id.drop_duplicates().values.tolist():
        pid = over_lapp[over_lapp['protein_id']==protein_id].sort_values(by=['Start location', 'Stop location'])
        start = pid['Start location'].values.tolist()
        stop = pid['Stop location'].values.tolist()
        pid = pid.protein_id.values.tolist()[0]
        c=0
        pids_o = []

        for s in start:

            if c>0 and c<(len(start)):

                if start[c] < stop[c-1]:
                    if pid not in pids_o:
                        pids_o.append(pid)
            c+=1

    print('OVERLAPPS of domain regions in these protein ids:')
    print(pids_o)

    # reshape the df
    interpro_df = reshape_interproscan_file(interpro_df)

    # get the domain ark
    domain_arc_list, cnt_domain_arc = get_domain_arc(interpro_df)

    # Sort the domain arks based on number of domains and domain location of shosen domain
    domain_arc_list.sort(key=dom_loc, reverse=True) # sort on chosen domain
    if sort_on_domain != []:
        domain_arc_list.sort(key=number_of_dom, reverse=True)

    schematic_df = create_schematic_table(domain_arc_list, cnt_domain_arc, interpro_df, interpro_original_df)
    # Save the df to a tsv file, this will be used to produce the domain arkitcture figure
    f = schematic_df.to_csv(out_file, header=None, index=False, sep='\t')
    return f


##### get input output and run the script

interpro_path = snakemake.input.interpro_path
out_file_path = snakemake.output.out_file_path
analysis = snakemake.params.analysis
merged_domains = snakemake.input.merged_domains
phobius = snakemake.input.phobius
pids = pd.read_csv(snakemake.input.pids, sep='\t', header=None)[0].values.tolist()
groups = pd.read_csv(str(snakemake.input.groups_protein_ids), sep='\t')
interproscan_header = ['protein_id','Sequence length','Analysis', 'Signature accession', 'Signature description', 'Start location', 'Stop location', 'Score','Date of run']
schematic = 'yes'
sort_on_domain = ['PF00069']
potential_domains = 100
global selected_domains
selected_domains = sort_on_domain


# run the script:
main(interpro_path, analysis, out_file_path, merged_domains, phobius, pids, groups, interproscan_header, schematic, potential_domains, sort_on_domain)



