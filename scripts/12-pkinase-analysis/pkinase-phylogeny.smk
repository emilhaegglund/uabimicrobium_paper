"""
Author: Emil Hagglund
Date: 240124

Workflow to perform phylogenetic analysis of the Pkinase domain inside the
PVC.
"""
rule all:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-clusters-90.faa",
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-clusters-80.faa",
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.bmge_blosum30.fasttree.nwk",
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.bmge_blosum30.fasttree.nwk",
        # Test other methtod
        #"../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-pfam-domains.txt",
        #"../../results/12-pkinase-analysis/pkinase-hmm/PF00069.28.hmm",
        #"../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-pfam-domains.stm"



rule extract_pkinase_seq:
    """
    Use inforrmation from Supplementary Table 5 to extract the protein
    sequences for the Pkinase domains from PVC genomes and write them to a
    fasta file.
    """
    input:
        suppl_table_5 = "../../tables/supplementary-table-5.tsv",
        proteoeme_dir = "../../results/01_download_data/genbank_proteomes/"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.faa"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python extract-pkinase-sequences.py {input.suppl_table_5} {input.proteoeme_dir} {output}"

rule cluster_pkinase_90:
    """
    Test clustering the Pkinase sequences using CD-HIT at 0.9 identity.
    Results: Clustered 2651 sequences to 2539 clusters.
    """
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.faa"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-clusters-90.faa"
    conda:
        "../envs/cd-hit.yaml"
    shell:
        "cd-hit -i {input} -o {output}"

rule cluster_pkinase_80:
    """
    Test clustering the Pkinase sequences using CD-HIT at 0.8 identity
    Results: Clustered 2651 sequences to 2396 clusters.
    """
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.faa"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-clusters-80.faa"
    conda:
        "../envs/cd-hit.yaml"
    shell:
        "cd-hit -c 0.8 -i {input} -o {output}"

rule extract_hmm_profile_from_pfam:
    """
    Extract the PF00069-hmm from Pfam database to use as a guide in the MSA
    """
    input:
        pfam = "/data/emil/pfam_v35.0/Pfam-A.hmm"
    output:
        "../../results/12-pkinase-analysis/pkinase-hmm/PF00069.28.hmm"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmfetch {input.pfam} PF00069.28 > {output}"

rule hmmalign_original:
    """
    Clutering the sequences did not significantly decrease the total number
    of Pkinases. Continue with all. Make use of the HMM of the Pkinase domain
    to align the seuqences.
    """
    input:
        pkinase_hmm = "../../results/12-pkinase-analysis/pkinase-hmm/PF00069.28.hmm",
        pkinase_seq = "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.faa"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.stm"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmalign -o {output} {input.pkinase_hmm} {input.pkinase_seq}"

rule convert_stm_to_faa:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.stm"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.aln"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "esl-reformat afa {input} > {output}"

rule bmge_trimming:
    """Trim the alignment using BMGE."""
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.bmge_blosum30.aln"
    conda:
        "../envs/bmge.yaml"
    shell:
        "bmge -i {input} -of {output} -t AA -m BLOSUM30 -h 0.7 -g 0.8"


rule fasttree:
    """First create an initial phylogeny using FastTree."""
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.bmge_blosum30.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align.bmge_blosum30.fasttree.nwk"
    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"

rule remove_long_branches:
    input:
        pkinase_seq = "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase.faa",
        seq_to_clean = "../../data/pkinase_clean.txt"
    output:
        pkinase_seq = "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-cleaned.faa"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python clean-fasta.py {input.pkinase_seq} {input.seq_to_clean} {output}"

rule hmmalign_original_cleaned:
    """
    Clutering the sequences did not significantly decrease the total number
    of Pkinases. Continue with all. Make use of the HMM of the Pkinase domain
    to align the seuqences.
    """
    input:
        pkinase_hmm = "../../results/12-pkinase-analysis/pkinase-hmm/PF00069.28.hmm",
        pkinase_seq = "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-cleaned.faa"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.stm"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmalign -o {output} {input.pkinase_hmm} {input.pkinase_seq}"

rule convert_stm_to_faa_cleaned:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.stm"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.aln"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "esl-reformat afa {input} > {output}"

rule convert_x_to_gaps:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.aln"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python convert-x-to-gaps.py {input} {output}"


rule bmge_trimming_cleaned:
    """Trim the alignment using BMGE."""
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.aln"
    conda:
        "../envs/bmge.yaml"
    shell:
        "bmge -i {input} -of {output} -t AA -m BLOSUM30 -h 0.7"

rule no_x:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.no-x.aln"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python convert-x-to-gaps.py {input} {output}"

rule gap_trim:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.no-x.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.no-x.gap.aln"
    conda:
        "../envs/bmge.yaml"
    shell:
        "bmge -i {input} -t AA -h 1 -w 1 -g 0.998 -of {output}"


rule remove_seq_with_many_gaps:
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.no-x.gap.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.no-x.gap.clean.aln"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python remove-seq-with-many-gaps.py  {input} {output}"

rule fasttree_cleaned:
    """First create an initial phylogeny using FastTree."""
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.no-x.bmge_blosum30.no-x.gap.clean.aln"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-hmm-align-cleaned.bmge_blosum30.fasttree.nwk"
    conda:
        "../envs/fasttree.yaml"
    shell:
        "fasttree -lg {input} > {output}"



# Test another method to perform the same analysis
rule extract_pkinase_proteins:
    """
    First extract the proteins for which we have detected a Pkinase domain in PVC.
    """
    input:
        suppl_table_5 = "../../tables/supplementary-table-5.tsv",
        proteoeme_dir = "../../results/01_download_data/genbank_proteomes/"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-proteins.faa"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python extract-pkinase-proteins.py {input.suppl_table_5} {input.proteoeme_dir} {output}"

rule pfam_scan:
    """
    Run Pfam-scan to search for the PF00069-domain.
    """
    input:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-proteins.faa"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-pfam-domains.txt"
    params:
        db="/data/emil/pfam_v35.0/"
    conda:
        "../envs/pfam_scan.yaml"
    threads:
        12
    shell:
        "pfam_scan.pl -fasta {input} -dir {params.db} -outfile {output} -cpu {threads}"

rule extract_pkinase_pfam_seq:
    """
    Using the hits from the Pfam-scan, extract the domains.
    """
    input:
        pfams = "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-pfam-domains.txt",
        proteome_dir = "../../results/01_download_data/genbank_proteomes/"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-pfam-domains.faa"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python extract-pkinase-pfam-sequences.py {input.pfams} {input.proteome_dir} {output}"

rule hmmalign:
    input:
        pkinase_hmm = "../../results/12-pkinase-analysis/pkinase-hmm/PF00069.28.hmm",
        pkinase_seq = "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-pfam-domains.faa"
    output:
        "../../results/12-pkinase-analysis/pkinase-phylogeny/pvc-pkinase-pfam-domains.stm"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmalign -o {output} {input.pkinase_hmm} {input.pkinase_seq}"
