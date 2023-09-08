# search for Pfam domain PF00069 (protein kinase domain) in all pvc species and 3 eukaryotes (+ ecoli with no STKs)
euk_accessions = ["GCA_000002985.3",
                  "GCA_000002945.2",
                  "GCA_000146045.2"]

pvc_accessions = []
with open("../../data/accession_numbers.csv") as f:
    for line in f:
        line = line.strip()
        accession = line.split(",")[0]
        if accession != "Assembly Accession":
            pvc_accessions.append(accession)

rule all:
    input:
        expand("../../results/12-pkinase-analysis/hmmsearch-euk/{euk_accession}-pkinase-reformat.tsv", euk_accession=euk_accessions),
        expand("../../results/12-pkinase-analysis/hmmsearch-pvc/{pvc_accession}-pkinase-reformat.tsv", pvc_accession=pvc_accessions)
        "../../snakemake-dags/01-hmm-pkinase-domain-search.png"

rule download_proteomes:
    """ Download proteomes for selected eukaryotes """
    output:
        temp("../../results/12-pkinase-analysis/euk-proteomes/{euk_accession}.zip")
    params:
        accession = "{euk_accession}"
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        "datasets download genome accession {params.accession} --include protein --filename {output}"

rule gunzip_proteomes:
    """ Unzip proteomes """
    input:
        '../../results/12-pkinase-analysis/euk-proteomes/{euk_accession}.zip'
    output:
        '../../results/12-pkinase-analysis/euk-proteomes/{euk_accession}.faa'
    params:
        accession = "{euk_accession}"
    shell:
       "unzip -p {input} ncbi_dataset/data/{params.accession}/protein.faa > {output}"

rule download_pfam_domain_profiles:
    """ Download the Pfam hmm-profiles """
    output:
        "../../results/12-pkinase-analysis/Pfam-hmms/Pfam-A.hmm.gz"
    shell:
        "wget -O {output} ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz"

rule gunzip_pfam_hmm_db:
    """ Unzip hmm-profiles """
    input:
        "../../results/12-pkinase-analysis/Pfam-hmms/Pfam-A.hmm.gz"
    output:
        "../../results/12-pkinase-analysis/Pfam-hmms/Pfam-A.hmm"
    shell:
        "gunzip -k -q {input} {output}"

rule get_pkinase_domain_profile:
    """  Get the hmm-profile for the Pkinase domain """
    input:
        "../../results/12-pkinase-analysis/Pfam-hmms/Pfam-A.hmm"
    output:
        "../../results/12-pkinase-analysis/pkinase-hmm/PF00069.hmm"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmfetch -o {output} {input} PF00069.28"

rule hmm_search_pvc:
    """ Search the PVC genomes using the Pkinase hmm-profile as query """
    input:
        hmm_model="../../results/12-pkinase-analysis/pkinase-hmm/PF00069.hmm",
        proteome="../../results/01_download_data/genbank_proteomes/{pvc_accession}.faa"
    output:
        "../../results/12-pkinase-analysis/hmmsearch-pvc/{pvc_accession}-pkinase.tsv"
    conda:
        "../envs/hmmer.yaml"
    log:
        "../../logs/hmmsearch/{pvc_accession}.log"
    shell:
        "hmmsearch --domtblout {output} {input.hmm_model} {input.proteome} >> {log}"

rule hmm_search_eukaryotes:
    """ Search the eukaryotic genomes using the Pkinase hmm-profile as query """
    input:
        hmm_model="../../results/12-pkinase-analysis/pkinase-hmm/PF00069.hmm",
        proteome="../../results/12-pkinase-analysis/euk-proteomes/{euk_accession}.faa"
    output:
        "../../results/12-pkinase-analysis/hmmsearch-euk/{euk_accession}-pkinase.tsv"
    conda:
        "../envs/hmmer.yaml"
    log:
        "../../logs/hmmsearch/{euk_accession}.log"
    shell:
        "hmmsearch --domtblout {output} {input.hmm_model} {input.proteome} >> {log}"

rule reformat_hmm_euk_output:
    input:
        "../../results/12-pkinase-analysis/hmmsearch-euk/{euk_accession}-pkinase.tsv"
    output:
        "../../results/12-pkinase-analysis/hmmsearch-euk/{euk_accession}-pkinase-reformat.tsv"
    shell:
        "sed 's/  */ /g' {input} > {output}"

rule reformat_hmm_pvc_output:
    input:
        "../../results/12-pkinase-analysis/hmmsearch-pvc/{pvc_accession}-pkinase.tsv"
    output:
        "../../results/12-pkinase-analysis/hmmsearch-pvc/{pvc_accession}-pkinase-reformat.tsv"
    shell:
        "sed 's/  */ /g' {input} > {output}"

rule create_dag_figure:
    output:
        "../../snakemake-dags/01-hmm-pkinase-domain-search.png"
    shell:
        'snakemake -s snakefiles/01_hmm_pkinase_domain_search.smk --forceall --rulegraph | dot -Tpng > {output}'
