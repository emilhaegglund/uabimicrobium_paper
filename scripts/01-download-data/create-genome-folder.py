"""
Script to move the genome data from the NCBI datasets folder to a genome specific folder
"""
import os
import sys
import shutil
import sys

ncbi_datasets_path = sys.argv[1]
output_file_path = sys.argv[2]
output_dir, accession = os.path.split(output_file_path)
accession = os.path.splitext(accession)[0]
for f in os.listdir(os.path.join(ncbi_datasets_path, "data", accession)):
    f_accession = "_".join(f.split("_")[:2])
    if f_accession == accession:
        genome_path = os.path.join(ncbi_datasets_path, "data", accession, f)
        shutil.copy(genome_path, output_file_path)
