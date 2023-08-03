#### PLACE THE DN, DS, DNDS VAALUES IN ONE TABLE
import pandas as pd


def open_data(data_file, column_names):
    data = pd.read_csv(data_file, sep = '\s+', header = None, names = column_names, na_values='')
    return data


def parse_data(my_data):
    
    parsed_table = []
    first_row=True
    
    for i in my_data.index:
        
        if first_row == True:
            species = int(my_data.iloc[i, 0])

        if  first_row==False and i < species: # since the index starts on 0 and the first row is the amount of species. If you have another amout of species chamge to number of species - 2 instead of 64
    
            first_gene_id = my_data.iloc[i, 0]
            second_gene_id = my_data.iloc[i+1, 0]
            pairewise_value = my_data.iloc[i+1, i]
            parsed_table.append([first_gene_id, second_gene_id, pairewise_value])
        
        first_row = False

    parsed_table = pd.DataFrame(parsed_table, columns = ['gene_1','gene_2','value'])
    print(parsed_table)
    return parsed_table



### Merge the dN, dS and dNdS file values into one table with the index replaced with the actual protein name

def merge_tables(dN_table, dS_table):

    merge_table = pd.merge(dN_table, dS_table, on=['gene_1','gene_2'])
    merge_table.rename(columns={'value_x':'dN', 'value_y':'dS'}, inplace=True)
    merge_table['dNdS'] = round(merge_table['dN']/merge_table['dS'], 4)
    return merge_table


### RUN FUNCTIONS ####

# 1. Create a list of column names. The CODEML files have empty values that are tricky to open with pandas. Solve this by making your column lable list
column_names = list(range(0,65)) ## OBS if you have more then 65 species change to that number

# 2. Open the 2ML.dN codeml result file and parse the file to a file format of our choosing namely [gene1, gene2, comparison value] 
dN = open_data(str(snakemake.input[0]), column_names)#str(snakemake.input[0])
dN_table = parse_data(dN)

# 3. Do the same for the 2ML.dS file
dS = open_data(str(snakemake.input[1]), column_names)
dS_table = parse_data(dS)

# 4. Merge the dN and dS parsed files such as [gene1, gene2, dN value, dS value, dN/dS] 
merged_table = merge_tables(dN_table, dS_table)

# 5. Save the output file with the labels [gene1, gene2, dN value, dS value, dNdS value]
merged_table.to_csv(str(snakemake.output[0]), index=False, sep='\t')

