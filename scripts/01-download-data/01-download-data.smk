"""
Download genome
"""

ACCESSIONS = []
with open("../../data/accession_numbers.csv") as f:
    for line in f:
        line = line.strip()
        accession = line.split(",")[0]
        if accession != "Assembly Accession":
            ACCESSIONS.append(accession)

rule all:
    input:
        expand("../../results/01_download_data/genbank_proteomes/{accession}.faa", accession=ACCESSIONS),
        expand("../../results/01_download_data/genbank_genomes/{accession}.fna", accession=ACCESSIONS),
        expand("../../results/01_download_data/genbank_genes/{accession}.cds.fna", accession=ACCESSIONS),
        "../../results/01_download_data/protein_mapping.tsv",
        "../../tables/genome_statistics.tsv"

rule get_accessions:
    input:
        "../../data/accession_numbers.csv"
    output:
        "../../results/01_download_data/accession_numbers.txt"
    shell:
        "awk -F',' '{{ print $1  }}' {input} | grep -v 'Assembly Accession' > {output}"

rule download_data:
    input:
        "../../results/01_download_data/accession_numbers.txt"
    output:
        "../../results/01_download_data/ncbi_dataset.zip"
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        "datasets download genome accession --inputfile {input} --include genome,protein,cds,seq-report --filename {output}"

rule unzip:
    input:
        "../../results/01_download_data/ncbi_dataset.zip"
    output:
        directory("../../results/01_download_data/ncbi_dataset/")
    params:
        outdir = "../../results/01_download_data/"
    shell:
        "unzip {input} -d {params.outdir}"

rule create_proteome_folder:
    input:
        "../../results/01_download_data/ncbi_dataset/"
    output:
        "../../results/01_download_data/genbank_proteomes/{accession}.faa"
    shell:
        "python create-proteome-folder.py {input} {output}"

rule create_genome_folder:
    input:
        "../../results/01_download_data/ncbi_dataset/"
    output:
        "../../results/01_download_data/genbank_genomes/{accession}.fna"
    shell:
        "python create-genome-folder.py {input} {output}"

rule create_gene_folder:
    input:
        "../../results/01_download_data/ncbi_dataset/"
    output:
        "../../results/01_download_data/genbank_genes/{accession}.cds.fna"
    shell:
        "python create-gene-folder.py {input} {output}"

rule protein_mapping:
    input:
        expand("../../results/01_download_data/genbank_proteomes/{accession}.faa", accession=ACCESSIONS)
    output:
        "../../results/01_download_data/protein_mapping.tsv"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python protein-mapping.py {input} {output}"

rule download_metadata:
    input:
        "../../results/01_download_data/accession_numbers.txt"
    output:
        "../../results/01_download_data/genome_metadata.tsv"
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        datasets summary genome accession --inputfile {input} --as-json-lines | \
        dataformat tsv genome --fields accession,organism-name,assminfo-level,assmstats-total-sequence-len,assmstats-gc-percent,annotinfo-featcount-gene-protein-coding > {output}
        """

rule compile_genome_stats:
    input:
        accessions = "../../data/accession_numbers.csv",
        metadata = "../../results/01_download_data/genome_metadata.tsv"
    output:
        "../../tables/genome_statistics.tsv"
    shell:
        "python compile-genome-statistics.py {input.accessions} {input.metadata} {output}"
