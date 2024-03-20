import pandas as pd

df = pd.read_csv('../../data/accession_numbers.csv')
accessions = df[df['group'].isin(['Brocadiaceae','Saltatorellus','Orphan'])]['assembly_accession'].to_list()

rule all:
	input:
		expand('../../results/04_interproscan/raw_results/{accession}.tsv', accession=accessions),
		expand('../../results/04_interproscan/reformated_results/{accession}.tsv', accession=accessions)



rule run_interproscan:
	input:
		'../../results/01_download_data/genbank_proteomes/{accession}.faa'
	output:
		'../../results/04_interproscan/raw_results/{accession}.tsv'
	shell:
		'interproscan.sh -i {input} -b {output}'

rule reformat_interproscan_results:
	input:
		'../../results/04_interproscan/raw_results/{accession}.tsv.tsv'
	output:
		'../../results/04_interproscan/reformated_results/{accession}.tsv'
	shell:
		'python reformat_results.py {input} {output}'
