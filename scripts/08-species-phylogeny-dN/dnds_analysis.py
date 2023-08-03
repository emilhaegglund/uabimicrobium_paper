########## calculate the mean, median and max of dN/dS values of selected Orthogroups
########## generate a statistics table, boxplot figure and Uab vs. (Saltatorellus/Anammox) dN value table

import statistics as s
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# paths
path_chosen_OG = '../../results/02_orthofinder/single_copy_orthogroups.tsv'
path_out_folder = '../../results/10_species_phylogeny_slow_evolving/'
accession_table = '../../data/accession_numbers.csv'

######## Functions

# return the name of the chosen orthogroups (OG)

def retrieve_OG(path_chosen_OG):
  data = pd.read_csv(path_chosen_OG, header = None, sep='\t')
  OG_ids = data.iloc[:,1]
  return OG_ids

###### Create a statistics table with mean, median, min and max dN/dS values

# open the codeml file into a dataframe
def open_codeml_file(OG, path_out_folder):
	try:
		codeml_data = pd.read_csv(str(path_out_folder) + 'CODEML/codeml_parsed/' + str(OG).strip() + '_dnds_table.tsv', sep='\t')
	except:
		exit('codeml data file could not be opened')
	return codeml_data

# compute statistics for a codeml data file
def calculate_statistics(codeml_data):
    # extract the dN and dS column from the Orthogroup file
    dN = codeml_data.iloc[:,2]
    print(dN)
    dS = codeml_data.iloc[:,3]
    print(dS)
    # save the mean, median and max value of the Orthogroup to the list data.
    data = [round(s.mean(dN), 5), round(s.median(dN),5), max(dN), min(dN), round(s.mean(dS),5), round(s.median(dS),5), max(dS), min(dS)]
    return data, dN

# Uses the two functions above to open codeml file and calculate statistics for all Orthogroups - return statistical table
def create_stat_table(OG_ids, path_out_folder):

	stat_table = []
	dN_list = []

	# calculate the statistics for all Orthogroups choosen and save to data frame
	for OG in OG_ids:
		codeml_file = open_codeml_file(OG, path_out_folder)
		data, dN = calculate_statistics(codeml_file)
		stat_table.append(data)
		dN_list.append(dN)

	# the column names for the computed statistics
	columns = ['dn_mean_value', 'dn_median', 'dn_max_value', 'dn_min_value', 'ds_mean_value', 'ds_median', 'ds_max_value', 'ds_min_value']

	#create the dataframe with the calculated statistical data and column names
	stat_data = pd.DataFrame(stat_table, index = OG_ids, columns = columns)

	return stat_data, dN_list


# Save the statistics table as a file
def save_statistics_table(stat_data, path_out_folder):

    f = stat_data.to_csv(str(path_out_folder) + 'dNdS_statistics.tsv', sep='\t')

    return f



#################### dN boxplot figure, sorted by dN ascending mean values

def mean_dN_boxplot(OG_ids, dN_list, path_out_folder):

	fig, ax = plt.subplots()
	ax.set_title('dN')

	# create dataframe with dN values and OG ids as columns
	df = pd.DataFrame(dN_list).T
	df.columns = OG_ids

	# sort the dataframe based on the dN mean values
	sorted_index = df.mean(axis = 0, skipna = True).sort_values().index
	df_sorted = df[sorted_index]

	# create the boxplot of sorted dN values of the Orthogroups
	sns.boxplot(data = df_sorted, showmeans = True, showfliers = False)

	# save figure with design settings (size and rotation of labels)
	fig.set_size_inches(15, 10)
	plt.xticks(rotation = 90)
	plt.savefig(str(path_out_folder) + 'dN_orthogroups.png')

	return fig

#################### check the dN for Uab compaired with saltatorella and anammox. To analyse how similar they are to one and another

def compair_Uab(OG_ids, path_out_folder, accession_numbers):
    # The anammox, saltatorellus and Uab amorphum ids
    accession_numbers = pd.read_csv(str(accession_numbers))

    anammox = list(accession_numbers[accession_numbers['Group'] == 'Brocadiaceae']['Assembly Accession'])
    saltatorella = list(accession_numbers[accession_numbers['Group'] == 'Saltatorellus']['Assembly Accession'])
    Uab = ['GCA_009002475.1']

    f_out = []

    for OG in OG_ids:
        with open(str(path_out_folder) + 'CODEML/codeml_parsed/' + str(OG).strip() + '_dnds_table.tsv', 'r') as f:
            for line in f:
                l = line.split('\t')
                genome_1 = str(l[0]).split('.')[0]
                genome_2 = str(l[1]).split('.')[0]

                if genome_1 in Uab:
                    if genome_2 in saltatorella:
                        f_out.append([str(OG), 'Saltatorella', str(genome_2).strip(), str(l[2]).strip()])
                    if genome_2 in anammox:
                        f_out.append([str(OG), 'Anammox', str(genome_2).strip(), str(l[2]).strip()])

                if genome_2 in Uab:
                    if genome_1 in saltatorella:
                        f_out.append([str(OG), 'Saltatorella', str(genome_1).strip(), str(l[2]).strip()])
                    if genome_1 in anammox:
                        f_out.append([str(OG), 'Anammox', str(genome_1).strip(), str(l[2]).strip()])
        f.close()
    return f_out

# save the dN comparision to a file
def save_comp_file(f_out, path_out_folder):

    df = pd.DataFrame(f_out, columns = ['Orthogroup','Clade','Taxon id','dN comparison with Uab'])
    f = df.to_csv(str(path_out_folder) + 'dN_Uab_placment.tsv', sep='\t')

    return f


### write the lowest dN mean OGs to a list

def lowest_and_fastest_mean_dN(amount, OG_ids, dN_list, path_out_folder):
    # create dataframe with dN values and OG ids as columns
    df = pd.DataFrame(dN_list).T
    df.columns = OG_ids

    # sort the dataframe based on the dN mean values
    sorted_index = df.mean(axis = 0, skipna = True).sort_values().index
    df_sorted = df[sorted_index]
    print(df_sorted)
    slowest_OG = df_sorted.columns[0:amount]
    fastest_OG = df_sorted.columns[-(amount+1):-1]
    print(fastest_OG)
    print(slowest_OG)

    return slowest_OG, fastest_OG

# Save the slowest evolving OG:s names to a txt file
def save_lowest_mean(slowest, path_out_folder):

    f = open(str(path_out_folder) + '/dN_low/orthogroups.txt', 'w')

    for line in slowest_OG:
        f.write(line + '\n')

    return f

# Save the fastest evolving OG:s names to a txt file
def save_highest_mean(fastest, path_out_folder):

    f = open(str(path_out_folder) + '/dN_high/orthogroups.txt', 'w')

    for line in fastest_OG:
        f.write(line + '\n')

    return f


############# Run the functions

# Store the chosen Orthogroup ids in OG_ids
OG_ids = retrieve_OG(path_chosen_OG)

# 1. create the statistics table of dn/dS information for each Orthogroup and save as file
stat_data, dN_list = create_stat_table(OG_ids, path_out_folder)
save_statistics_table(stat_data, path_out_folder)

# 2. create a boxplot figure sorted by the mean dN values for each Orthogroup
mean_dN_boxplot(OG_ids, dN_list, path_out_folder)

# 3. Compair the Saltatorellus and the Annamox bacterial genomes' dN values with Uab amorphum and save to file
f_out = compair_Uab(OG_ids, path_out_folder, accession_table)
save_comp_file(f_out, path_out_folder)

# 4. extract the 20 slowest and fastest evolving Orthogroups and save the OG:s to a file
slowest_OG, fastest_OG = lowest_and_fastest_mean_dN(25, OG_ids, dN_list, path_out_folder)
save_lowest_mean(slowest_OG, path_out_folder)
save_highest_mean(fastest_OG, path_out_folder)
