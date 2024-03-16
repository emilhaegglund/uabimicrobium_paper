import pandas as pd
def read_interpro(interpro_path):
    col_names = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession', 'Signature description', 'Start location', 'Stop location', 'Score', 'Status',
        'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
        'optional: GO terms', 'optional: Pathway annotation']
    interproscan_df = pd.read_csv(interpro_path, sep='\t', header=None, names=col_names)
    return interproscan_df

def create_comparison_table(domain_table,interpro_filtered,outfile):
    hmm_active_df = pd.read_csv(domain_table, sep='\t')
    hmm_active_pid = hmm_active_df['protein_id'].drop_duplicates().values.tolist()

    pid_interpro_df = read_interpro(interpro_filtered)
    pid_interpro = pid_interpro_df[pid_interpro_df['Signature accession']=='PF00069']['protein_id'].drop_duplicates().values.tolist()

    p_id = []
    for pid in hmm_active_pid:
        if pid not in pid_interpro:
            p_id.append(pid)

    comp_df = hmm_active_df[hmm_active_df['protein_id'].isin(p_id)]
    comp_df.to_csv(outfile, index=False, sep='\t')
    return comp_df



domain_table = str(snakemake.input.domain_table)
interpro_filtered = str(snakemake.input.interpro_filtered)
outfile = str(snakemake.output.outfile)

df = create_comparison_table(domain_table,interpro_filtered,outfile)