# Compare all proteins to the COG database ausing blast and assign functional annotation

rule all:
	input:
		'../../results/05_cog/all_proteins.COG'



rule merge_proteomes:
	output:
		'../../results/05_cog/all_proteins.fasta'
	script:
		'merge_proteomes.py'

rule run_cog:
	input:
		'../../results/05_cog/all_proteins.fasta',
		'/home/dtamarit/cog170403/prot2003-2014_formatted.fa'
	output:
		'../../results/05_cog/all_proteins_vs_cog_db.blastp'
	threads:
		40

	#envs:
	shell:
		'blastp -query {input[0]} -db {input[1]} -evalue 1e-5 -num_threads {threads} > {output}'

rule assign_cog_categories:
	input:
		'../../results/05_cog/all_proteins_vs_cog_db.blastp'
	output:
		'../../results/05_cog/all_proteins.COG.raw'
	shell:
		'perl COG_assignment.pl {input} > {output}'

rule parse_cog_results:
	input:
		'../../results/05_cog/all_proteins.COG.raw'
	output:
		'../../results/05_cog/all_proteins.COG'
	shell:
		'perl parse_COG_assignment.pl {input} > {output}'
