"""
Script to move the proteome data from the NCBI datasets folder to a proteome specific folder
"""
import os
import sys
import shutil
import sys

ncbi_datasets_path = sys.argv[1]
output_file_path = sys.argv[2]
output_dir, accession = os.path.split(output_file_path)
print(accession)
accession = os.path.splitext(accession)[0]
print(accession)

proteome_path = os.path.join(ncbi_datasets_path, "data", accession, "protein.faa")
shutil.copy(proteome_path, output_file_path)





