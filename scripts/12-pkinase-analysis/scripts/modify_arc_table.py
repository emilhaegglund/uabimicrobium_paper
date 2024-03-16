import pandas as pd



arc_table_path = snakemake.input.arc_table_path
out_file_path = snakemake.output.out_file_path

df = pd.read_csv(arc_table_path, sep='\t')

print(df)