#!/usr/bin/env python
"""
Filter blast results based on the alignment length. Used in the alignment step of the clustering of protein families.

Run:
./alignment_length_filter.py blast_in.tsv filtered_blast_out.tsv

Described in https://doi.org/10.1093/molbev/msz287
"""
import pandas as pd
import sys
import numpy as np


blast_format_six_headers = [
    "qseqid",
    "sseqid",
    "qlen",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]
blast_df = pd.read_csv(
    sys.argv[1],
    sep="\t",
    header=None,
    names=blast_format_six_headers,
    compression="gzip",
)

blast_df["percentage"] = (5 * np.log(1600 / blast_df["qlen"]) / np.log(2)) + 45
mask = blast_df.percentage < 45
blast_df.loc[mask, "percentage"] = 45

mask = blast_df.percentage > 70
blast_df.loc[mask, "percentage"] = 70

blast_df["filter_length"] = (blast_df["qlen"] * blast_df["percentage"]) / 100
blast_df["qalignlen"] = abs(blast_df["qend"] - blast_df["qstart"])
blast_df = blast_df[blast_df["qalignlen"] >= blast_df["filter_length"]]
blast_df[
    [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
].to_csv(sys.argv[2], sep="\t", index=False, header=False)
