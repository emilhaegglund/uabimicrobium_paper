import pandas as pd

def read_table(table_path):
    df = pd.read_csv(table_path, sep='\t', header=None)
    df = df[[0,2,5,6,7]]
    df.columns = ['Architecture','Domains','# of proteins','protein ids','Domain annotations']
    return df


def join_rows(df, grouped_by, joined_column):
    df = df.groupby(grouped_by).agg({joined_column: lambda x: ';'.join(x)})
    df = df.reset_index()
    return df

def reformat_table(df, df_anno_domains, df_anno_protein, df_group, df_og):
    df = df[['Architecture','# of proteins','protein ids']]
    df = df.iloc[::-1]

    df = df.join(df_anno_protein.set_index('Architecture'), on='Architecture')
    df = df.join(df_anno_domains.set_index('Architecture'), on='Architecture')
    df = df.join(df_og.set_index('Architecture'), on='Architecture')
    df = df.join(df_group.set_index('Architecture'), on='Architecture')
    df = df.drop_duplicates()

    return df

def extract_group_ids_per_arc(df, table_path):
    df_g = pd.read_csv(table_path, sep='\t')
    df_g = df_g[['group','pkinase protein ids']]
    df_g = df_g.groupby('group').agg({'pkinase protein ids': lambda x: ';'.join(x)})
    df_g = df_g.reset_index()

    groups_to_add = []
    for group in df_g.group:
        pid_group = df_g[df_g['group']==group]['pkinase protein ids'].values.tolist()[0].split(';')

        for Arc in df.Architecture.drop_duplicates().values.tolist():
            pid_arc = df[df['Architecture']==Arc]['protein ids'].values.tolist()[0].split(';')

            for pid in pid_arc:
                if pid in pid_group:
                    groups_to_add.append([Arc, group])


    df_group_to_add = pd.DataFrame(groups_to_add)

    df_group_to_add = df_group_to_add.groupby(0).agg({1: lambda x: ';'.join(x)})
    df_group_to_add = df_group_to_add.reset_index()
    df_group_to_add.columns = ['Architecture','Group']

    return df_group_to_add

def read_og_file(og_file_path, clustering_software):
    og_df = pd.read_csv(og_file_path, sep='\s', header=None)
    og_df = og_df.melt(id_vars=0).reset_index(drop=True)
    og_df = og_df.dropna() # remove lines with missing values after the melt command
    og = og_df[[0, 'value']] # extract the orthogroup and protein columns
    og.columns = [ str(clustering_software) + '_orthogroup', 'sseqid'] # rename the orthogroup and protein columns (sseqid will be used as key for mergin to blast results)
    return og

def add_protein_anno_og(df, anno_path, og_df_orthomcl):

    df_protein_anno = pd.read_csv(anno_path, sep='\t')
    df_protein_anno = df_protein_anno.fillna('Missing')
    df_protein_anno = df_protein_anno[['annotation', 'pkinase protein ids']]
    data_to_add = []

    for row_index in df_protein_anno.index:
        row = df_protein_anno.iloc[row_index,:]

        anno = row[0].split(';')
        pids = str(row[1]).split(';')

        for i in range(len(pids)-1):
            data_to_add.append([anno[i], pids[i]])

    df_protein_anno = pd.DataFrame(data_to_add)
    df_protein_anno.columns = ['Protein annotation', 'Protein id']

    arc_og_data = []
    arc_anno_data = []
    for Arc in df.Architecture.drop_duplicates().values.tolist():
        pid_arc = df[df['Architecture']==Arc]['protein ids'].values.tolist()[0].split(';')

        for pid in pid_arc:
            og = og_df_orthomcl[og_df_orthomcl['sseqid']==str(pid)]['orthomcl_orthogroup'].values.tolist()
            if og == []:
                og = 'Missing'
            else:
                og = og[0]
            anno_prot = df_protein_anno[df_protein_anno['Protein id']==str(pid)]['Protein annotation'].values.tolist()

            if anno_prot == []:
                anno_prot = 'Missing'
            else:
                anno_prot = anno_prot[0]

            arc_og_data.append([Arc, og])
            arc_anno_data.append([Arc, anno_prot])


    df_og_to_add = pd.DataFrame(arc_og_data)
    df_og_to_add = df_og_to_add.groupby(0).agg({1: lambda x: ';'.join(x)})
    df_og_to_add = df_og_to_add.reset_index()
    df_og_to_add.columns = ['Architecture','Orthomcl OG']

    df_anno_to_add = pd.DataFrame(arc_anno_data)
    print(df_anno_to_add)
    df_anno_to_add = df_anno_to_add.groupby(0).agg({1: lambda x: ';'.join(x)})
    df_anno_to_add = df_anno_to_add.reset_index()
    df_anno_to_add.columns = ['Architecture','Protein annotation']

    return df_anno_to_add, df_og_to_add

def main(table_path, kinase_path, og_path, out_file):
    df = read_table(table_path)
    og_df_orthomcl = read_og_file(og_path, 'orthomcl')

    df_domains = join_rows(df, 'Architecture', 'Domains')
    df_anno_domain = join_rows(df,'Architecture','Domain annotations')
    df_anno_protein, df_og = add_protein_anno_og(df, kinase_path, og_df_orthomcl)
    df_group = extract_group_ids_per_arc(df, kinase_path)

    df = reformat_table(df, df_anno_domain, df_anno_protein, df_group, df_og)
    f = df.to_csv(out_file, sep='\t', index=False)
    return f

## input/output
table_path = snakemake.input.table_path
kinase_path = snakemake.input.kinase_path
og_path = snakemake.input.og_path
out_file = snakemake.output.out_file

## run functions
f = main(table_path, kinase_path, og_path, out_file)
