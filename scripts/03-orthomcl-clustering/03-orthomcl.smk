import os
import random
import string

def get_names():
    """
    Generate a set of random names to use for the fasta in the orthoMCL pipeline
    """
    random.seed(0)
    names = []
    cwd = os.getcwd()
    num_files = len([f for f in os.listdir('../../results/01_download_data/genbank_proteomes/') if os.path.splitext(f)[1] == '.faa'])
    while len(names) < num_files:
        name = "".join(random.choice(string.ascii_uppercase) for _ in range(4))
        if name not in names:
            names.append(name)
    return names

NAMES = get_names()

rule all:
    input:
        expand("../../results/03_orthomcl_clustering/orthomcl_faa/{name}.fasta", name=NAMES),
        "../../results/03_orthomcl_clustering/orthomcl_run/all_proteins.fasta",
        "../../results/03_orthomcl_clustering/orthomcl_run/all_proteins.dmnd",
        "../../results/03_orthomcl_clustering/orthomcl_run/diamond_all_vs_all.tab.gz",
        "../../results/03_orthomcl_clustering/orthomcl_run/diamond_all_vs_all.filter.tab",
        "../../results/03_orthomcl_clustering/install_schema.log",
        "../../results/03_orthomcl_clustering/orthomcl_run/Orthogroups.translated.txt",
        "../../results/03_orthomcl_clustering/orthomcl_results/Orthogroups.txt"

rule setupDB:
    """
    Setup the MySQL database.
    """
    input:
        orthomcl_config="orthomcl-config-template.txt"
    output:
        "../../results/03_orthomcl_clustering/install_schema.log"
    conda:
        "../envs/orthomcl.yaml"
    shell:
        """
        mysql -u orthomcl_user -p -e "create database orthomcl_2023";
        orthomclInstallSchema {input.orthomcl_config} {output}
        """

rule compliantFasta:
    """
    Code to replace the orthomclAdjustFasta step because the
    format of input didn't fit well with the snakemake script.
    The code snippet does the same as the orthomclAdjustFasta
    script, changes the name of the file and the sequence headers.
    """
    output:
        seqs=expand("../../results/03_orthomcl_clustering/orthomcl_faa/{name}.fasta", name=NAMES),
        species_dict="../../results/03_orthomcl_clustering/orthomcl_run/species_dict.tsv",
        proteome_dict="../../results/03_orthomcl_clustering/orthomcl_run/proteome_dict.tsv"
    run:
        import os
        from Bio import SeqIO

        if not os.path.exists("../../results/03_orthomcl_clustering/orthomcl_faa"):
            os.mkdir("../../results/03_orthomcl_clustering/orthomcl_faa")
        species_dict = open(output.species_dict, 'w')
        proteome_dict = open(output.proteome_dict, 'w')
        for i, proteome in enumerate([f for f in os.listdir('../../results/01_download_data/genbank_proteomes/') if os.path.splitext(f)[1] == '.faa']):
            new_headers = []
            proteome_path = os.path.join('../../results/01_download_data/genbank_proteomes', proteome)
            records = SeqIO.parse(proteome_path, 'fasta')
            pre_id = sorted(NAMES)[i]
            species_dict.write(pre_id + '\t' + proteome + '\n')
            for record in records:
                old_part = record.id.split(' ')[0]
                counter = old_part.split('_')[-1]
                new_record_id = pre_id + '|' + pre_id + '_' + counter
                proteome_dict.write(new_record_id + '\t' + record.id + '\t' + record.description + '\n' )
                record.description = ''
                record.id = new_record_id
                new_headers.append(record)
            print(pre_id)
            SeqIO.write(new_headers, '../../results/03_orthomcl_clustering/orthomcl_faa/' + pre_id + '.fasta', 'fasta')
        species_dict.close()
        proteome_dict.close()

rule combineFasta:
    input:
        compliant_fasta=expand("../../results/03_orthomcl_clustering/orthomcl_faa/{name}.fasta", name=NAMES)
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/all_proteins.fasta",
    shell:
        "cat {input} > {output}"

rule diamondDB:
    input:
        "../../results/03_orthomcl_clustering/orthomcl_run/all_proteins.fasta"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/all_proteins.dmnd"
    conda:
        "../envs/diamond.yaml"
    shell:
        "diamond makedb --in {input} -p {threads} -d ../../results/03_orthomcl_clustering/orthomcl_run/all_proteins"

rule diamond:
    input:
        query="../../results/03_orthomcl_clustering/orthomcl_faa/{name}.fasta",
        db="../../results/03_orthomcl_clustering/orthomcl_run/all_proteins.dmnd"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/diamond_result/diamond_{name}.tab.gz"
    conda:
        "../envs/diamond.yaml"
    threads:
        4
    shell:
        """
        diamond blastp -q {input.query} --more-sensitive --db {input.db} --threads {threads} \
        --max-target-seqs 0 --outfmt 6 qseqid sseqid qlen pident length mismatch gapopen \
        qstart qend sstart send evalue bitscore --out {output} --compress 1
        """

rule concatDiamondResults:
    input:
        expand("../../results/03_orthomcl_clustering/orthomcl_run/diamond_result/diamond_{name}.tab.gz", name=NAMES)
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/diamond_all_vs_all.tab.gz"
    shell:
        "cat {input} > {output}"

rule filter_blast:
    input:
        "../../results/03_orthomcl_clustering/orthomcl_run/diamond_all_vs_all.tab.gz"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/diamond_all_vs_all.filter.tab"
    shell:
        "python alignment-length-filter.py {input} {output}"

rule parseBlast:
    input:
       blast_result="../../results/03_orthomcl_clustering/orthomcl_run/diamond_all_vs_all.filter.tab",
       fasta_files = expand("../../results/03_orthomcl_clustering/orthomcl_faa/{name}.fasta", name=NAMES)
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/similar_sequences.txt"
    conda:
        "../envs/orthomcl.yaml"
    shell:
        "orthomclBlastParser {input.blast_result} ../../results/03_orthomcl_clustering/orthomcl_faa/ >> {output}"

rule loadBlast:
    input:
        db_conf="orthomcl-config-template.txt",
        similar_sequences="../../results/03_orthomcl_clustering/orthomcl_run/similar_sequences.txt"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/loadedBlast"
    conda:
        "../envs/orthomcl.yaml"
    shell:
        """
        orthomclLoadBlast {input.db_conf} {input.similar_sequences};
        touch {output}
        """

rule pairs:
    input:
        finished_loadedBlast="../../results/03_orthomcl_clustering/orthomcl_run/loadedBlast",
        db_conf="orthomcl-config-template.txt"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/orthomcl_pairs.log"
    conda:
        "../envs/orthomcl.yaml"
    shell:
        "orthomclPairs {input.db_conf} {output} cleanup=no"

rule dumpPairsFiles:
    input:
        finished_pairs="../../results/03_orthomcl_clustering/orthomcl_run/orthomcl_pairs.log",
        db_conf="orthomcl-config-template.txt"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/mclInput"
    conda:
        "../envs/orthomcl.yaml"
    shell:
        """
        orthomclDumpPairsFiles {input.db_conf};
        mv mclInput {output}
        """

rule mcl:
    input:
        "../../results/03_orthomcl_clustering/orthomcl_run/mclInput"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/mclOutput"
    conda:
        "../envs/orthomcl.yaml"
    shell:
        "mcl {input} --abc -I 1.5 -o {output}"

rule groups:
    input:
        "../../results/03_orthomcl_clustering/orthomcl_run/mclOutput"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_run/Orthogroups.txt"
    conda:
        "../envs/orthomcl.yaml"
    shell:
        "orthomclMclToGroups cluster 1000 < {input} > {output}"

rule translate_back:
    input:
        orthogroups="../../results/03_orthomcl_clustering/orthomcl_run/Orthogroups.txt",
        proteome_dict="../../results/03_orthomcl_clustering/orthomcl_run/proteome_dict.tsv"
    output:
        orthogroups="../../results/03_orthomcl_clustering/orthomcl_run/Orthogroups.translated.txt"
    shell:
        "python get-original-protein-headers.py {input.orthogroups} {input.proteome_dict} {output}"

rule compile_orthomcl_resutls:
    input:
        orthogroups="../../results/03_orthomcl_clustering/orthomcl_run/Orthogroups.txt",
        species_dict="../../results/03_orthomcl_clustering/orthomcl_run/species_dict.tsv",
        proteome_dict="../../results/03_orthomcl_clustering/orthomcl_run/proteome_dict.tsv",
        protein_mapping="../../results/01_download_data/protein_mapping.tsv"
    output:
        "../../results/03_orthomcl_clustering/orthomcl_results/Orthogroups.txt",
        "../../results/03_orthomcl_clustering/orthomcl_results/Orthogroups.GeneCount.tsv"
    shell:
        """
        python orthomcl-compile-results.py --orthomcl {input.orthogroups} \
        --species_dictionary {input.species_dict} \
        --proteome_dictionary {input.proteome_dict} \
        --proteome ../../results/03_orthomcl_clustering/orthomcl_faa/ \
        --protein_mapping {input.protein_mapping} \
        --fraction 0.95 \
        --output ../../results/03_orthomcl_clustering/orthomcl_results/
        """
