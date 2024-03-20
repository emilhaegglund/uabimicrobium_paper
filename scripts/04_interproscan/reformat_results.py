"""
Written by Emil Hägglund

Script to reformat the output from InterProScan.
Output is a .tsv file where each line represents
a protein and the different annotations from the
different softwares. If a protein contains more
than one annotation from the same software, they
are separated with a ;.
"""
import sys
import pandas as pd


interpro_dict = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.strip('\n')
        line = line.split('\t')
        protein_id = line[0]
        software = line[3]
        name = line[5]
        if protein_id in interpro_dict.keys():
            interpro_dict[protein_id][software] += name + ';'
        else:
            interpro_dict[protein_id] = {'Pfam':'', 'TIGRFAM':'', 'SUPERFAMILY':'', 'Gene3D':'','CDD':'', 'SFLD':'', 'SMART':'', 'Hamap':'', 'ProSiteProfiles':'', 'MobiDBLite':'', 'Coils':'', 'PRINTS':'', 'ProSitePatterns':'', 'PIRSF':'', 'InterPro':'', 'InterPro_text':'', 'Domain start':'', 'Domain stop':''}
            interpro_dict[protein_id][software] += name + ';'
        if len(line) > 11:
            interpro_id = line[11]
            interpro_text = line[12]
            domain_start = line[6]
            domain_stop = line[7]
            interpro_dict[protein_id]['InterPro'] += interpro_id + ';'
            interpro_dict[protein_id]['InterPro_text'] += interpro_text + ';'
            interpro_dict[protein_id]['Domain start'] += domain_start + ';'
            interpro_dict[protein_id]['Domain stop'] += domain_stop + ';'

df = pd.DataFrame.from_dict(interpro_dict)
df = df.T
df['protein_id'] = df.index
df.reset_index(drop=True, inplace=True)
protein_id_col = df.pop('protein_id')
df.insert(0, 'protein_id', protein_id_col)

df.to_csv(sys.argv[2], index=False, sep='\t')
