import os
import sys
from Bio import SeqIO
import pandas as pd


# Read supplementary table which lists the positions of the Pkinase proteins
# and their positions in the protein.
col_names = ["protein_id", "domain_id",	"domain_from", "domain_to",	"protein_from",	"protein_to",	"catalytic_site_lys", "catalytic_site_asp_1",	"catalytic_site_asp_2",	"RD_site",	"domain_completness",	"group"]
df = pd.read_csv(sys.argv[1], sep="\t", names=col_names)
pkinase_list = df["protein_id"].unique().tolist()  # Extract pkinase proteins

# Loop through the proteomes and extract the sequences for the Pkinase domain
f_out = open(sys.argv[3], "w")
for f in os.listdir(sys.argv[2]):
    f_path = os.path.join(sys.argv[2], f)
    print(f)
    for record in SeqIO.parse(f_path, "fasta"):
        if record.id in pkinase_list:
            start_positions = df[df["protein_id"] == record.id]["protein_from"].tolist()
            stop_positions = df[df["protein_id"] == record.id]["protein_to"].tolist()
            # Consider that a single protein can have multiple Pkinase domains
            for j in range(len(start_positions)):
                    start_pos = int(start_positions[j])
                    stop_pos = int(stop_positions[j])
                    f_out.write('>' + record.id + "_" + str(j + 1) + "\n")
                    f_out.write(str(record.seq[start_pos:stop_pos]) + "\n")
f_out.close()