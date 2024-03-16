import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

domain_table = snakemake.input.domain_table
domain_list = snakemake.params.domain_list
accession = snakemake.params.accession
out_folder = snakemake.params.out_folder


def main(domain_table, domain_list, accession, out_folder):
    
    for domain in domain_list:
        df = pd.read_csv(domain_table, sep='\t')

        # get uab domain (tpr and wd40) counts to make #domains vs # of proteins figures
        df_domain_uab = df[df['accession']==accession][['specie', 'protein_id', str(domain)]]
        df_domain_uab = df_domain_uab[df_domain_uab[str(domain)]!=0][['specie',str(domain)]].value_counts()
        df_domain_uab = df_domain_uab.to_frame().reset_index()
        df_domain_uab.columns = ['specie', '# ' + str(domain), '# proteins']

        # save the figure
        ax = sns.barplot(x='# ' + str(domain), y="# proteins", hue="specie", data=df_domain_uab, errwidth=0)
        f = ax.get_figure()
        print('../results/uab_barplot_' + str(domain).replace(' ', '_') +'.png')
        f.savefig('../results/uab_barplot_' + str(domain).replace(' ', '_') +'.png')
        plt.close(f)

        # get pvc (non-uab) domain (tpr and wd40) counts to make #domains vs # of proteins figures
        df_domain_pvc_non_uab = df[~df['accession'].isin([accession, 'GCA_000002985.3', 'GCA_000002945.2', 'GCA_000146045.2'])][['specie', 'protein_id', str(domain), 'group']]

        # create figures for each group
        for g in df['group'].drop_duplicates().values.tolist():
            d = df[df[str(domain)]!=0]
            d = d[d['group'] == g][['specie', str(domain)]]

            if not d.empty:
                d = d.value_counts()
                d = d.to_frame().reset_index()

                d.columns = ['specie', '# ' + str(domain), '# proteins']
                ax = sns.barplot(x='# ' + str(domain), y="# proteins", hue="specie", data=d, errwidth=0)
                f = ax.get_figure()
                print('../results/' + str(g) + '_barplot_' + str(domain).replace(' ', '_') +'.png')
                f.savefig('../results/' + str(g) + '_barplot_' + str(domain).replace(' ', '_') +'.png')
                plt.close(f)  
    return f

f = main(domain_table, domain_list, accession, out_folder)
#ax = sns.stripplot(x="sex", y="total_bill", hue="day", data=tips)