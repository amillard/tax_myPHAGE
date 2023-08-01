#!/usr/bin/env python3
import pandas as pd
import subprocess
import io
from icecream import ic
import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import shutil
import sys
from argparse import ArgumentParser
import gzip


# Define color codes
class TerminalColors:
    # Regular text colors
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'

    # Text colors with a background
    BLACK_BG = '\033[40m'
    RED_BG = '\033[41m'
    GREEN_BG = '\033[42m'
    YELLOW_BG = '\033[43m'
    BLUE_BG = '\033[44m'
    MAGENTA_BG = '\033[45m'
    CYAN_BG = '\033[46m'
    WHITE_BG = '\033[47m'

    # Special text styles
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# this is the location of where the script and the databases are (instead of current_directory which is the users current directory)
HOME = os.path.dirname(__file__)


#check the files are present
VMR_path = os.path.join(HOME, 'VMR.xls')

blastdb_path = os.path.join(HOME, 'Bacteriophage_genomes.fasta')
ICTV_path = os.path.join(HOME, 'ICTV.msh')

if os.path.exists(VMR_path):
    print (TerminalColors.BLUE +f"Found {VMR_path} as expected")
else:
    print(TerminalColors.RED +f'File {VMR_path} does not exist')

if os.path.exists(ICTV_path):
    print(TerminalColors.BLUE +f" Found {ICTV_path} as expected")
else:
    print(TerminalColors.RED+f'File {ICTV_path} does not exist, was it downloaded correctly?')

#search for viridic before it is run
path_viridic_bash = os.path.join(HOME,'viridic_v1.1','viridic.bash')

if os.path.exists(path_viridic_bash):
    print (TerminalColors.BLUE +f"Found {path_viridic_bash} as expected")
else:
    print(TerminalColors.RED + f"""File {path_viridic_bash} does not exist 
 Expecting to find VIRIDIC installed. Need to install it or change the path if installed elsewhere""")


def is_program_installed_unix(program_name):
    try:
        subprocess.check_output(f'which {program_name}', shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

# check programs are installed
program_name = "blastdbcmd"
if is_program_installed_unix(program_name):
    ic(f"{program_name} is installed will proceed ")
else:
    print(TerminalColors.RED +f"{program_name} is not installed.")
    exit()

# check programs are installed
program_name = "mash"
if is_program_installed_unix(program_name):
    ic(f"{program_name} is installed will proceed ")
else:
    print(TerminalColors.RED +f"{program_name} is not installed.")
    exit()



#check if blastDB is present
if os.path.exists(blastdb_path):
    print (TerminalColors.BLUE +f"Found {blastdb_path} as expected")
else:
    print(TerminalColors.RED +f"File {blastdb_path} does not exist will create database now  ")
    print(TerminalColors.RED + "Will download the database now and create database")
    url = "https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz"
    file_download = "Bacteriophage_genomes.fasta.gz"
    extracted_file_name = "Bacteriophage_genomes.fasta"
    download_command = f"curl -o {file_download} {url}"
    try:
        subprocess.run(download_command, shell=True, check=True)
        print(f"{url} downloaded successfully!")
    except subprocess.CalledProcessError as e:
            print(f"An error occurred while downloading {url}: {e}")
    #Gunzip the file
    gunzip_command = f"gunzip {file_download}"
    try:
        subprocess.run(gunzip_command, shell=True, check=True)
        print("File gunzipped successfully!")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while gunzipping the file: {e}")

    # Run makeblastdb
    makeblastdb_command = f"makeblastdb -in {extracted_file_name} -parse_seqids -dbtype nucl"
    try:
        subprocess.run(makeblastdb_command, shell=True, check=True)
        print("makeblastdb command executed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing makeblastdb: {e}")


# Create the argument parser
usage = "%prog [options] file (or - for stdin)"
description='Read input.fasta file'
parser = ArgumentParser(usage, description=description)
parser.add_argument("-v", "--verbose", action="store_true", default = 0)
parser.add_argument("-t", "--threads", dest='threads', type=str, default= "8" ,
                    help= "Maximum number of threads that will be used")
parser.add_argument('-i', '--input', dest='in_fasta', type=str, help='Path to an input fasta file')
parser.add_argument("-p", "--prefix",  type=str, default ="", dest='prefix',
                    help='will add the prefix to results and summary files that will store results of MASH and comparision to the VMR Data produced by'
                         'ICTV combines both sets of this data into a single csv file. '
                         'Use this flag if you want to run multiple times and keep the results files ')


args, nargs = parser.parse_known_args()
verbose = args.verbose

#turn on ICECREAM reporting
if not verbose: ic.disable()


#Defined and set some parameters
query = os.path.join(HOME,'query.fasta')
fasta_file = args.in_fasta
threads = args.threads
ic(f"Number of set threads {threads}")
#create results folder
results_path = os.path.join(HOME,"results")

#path to the combined df containing mash and VMR data
out_csv_of_taxonomy = args.prefix+"Output_of_taxonomy.csv"
taxa_csv_output_path = os.path.join(results_path, out_csv_of_taxonomy)

#path the final results summary file
summary_results = args.prefix+"Summary_file.txt"
summary_output_path = os.path.join(results_path, summary_results)


if not os.path.exists(results_path):
    os.makedirs(results_path)
    print(TerminalColors.BLUE +f"Directory '{results_path}' created to store results")
else:
    print(TerminalColors.BLUE +f"\n\nWarning: Directory '{results_path}'already exists. All results will be overwritten.")


#fasta file to store known taxa
known_taxa_path = os.path.join(results_path, 'known_taxa.fa')
#store files for VIRIDIC run
viridic_in_path = os.path.join(results_path, 'viridic_in.fa')
#viridic path
viridic_path = os.path.join("")


#write the input file to a new file with a header called "taxmyPhage" which makes it easier to find in data later on

handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file)
entries = SimpleFastaParser(handle)

with open(fasta_file, "r") as input_file:
    records = list(SeqIO.parse(args.in_fasta, "fasta"))
    # Check the number of sequences in the file
    num_sequences = len(records)
    if num_sequences == 0:
        print(TerminalColors.RED+"Error: The FASTA file is empty.")
        exit()
    elif num_sequences == 1:
                # Open output FASTA file
        with open(query, "w") as output_file:
            #for name, seq in entries:
            for record in records:
                with open(summary_output_path, 'a') as file:
                    file.write(f"Query sequence header was:{record.description}\n")
                print(f">taxmyPhage\n{record.seq}", file=output_file)
    else:
        print(TerminalColors.RED +f"\nError: The {args.in_fasta} FASTA file contains {num_sequences} sequences."\
              f" Only one sequence can be classified at a time ")

#Read the viral master species record into a DataFrame
taxa_df = pd.read_excel(VMR_path,sheet_name=0)

#Print the DataFrame and rename a column
ic(f"taxa_df")

taxa_df = taxa_df.rename(columns={'Virus GENBANK accession': 'Genbank'})
taxa_df['Genbank'].fillna('', inplace=True)
# Get the column headings as a list

headings = list(taxa_df.columns)
#create a dictionary of Accessions linking to Genus
accession_genus_dict = taxa_df.set_index('Genbank')['Genus'].to_dict()

#Print the column headings one per line- error checking
for heading in headings:
    ic(heading)

#run mash to get top hit and read into a pandas dataframe

mash_output = subprocess.check_output(['mash', 'dist', '-d' , '0.2',  '-p', threads, ICTV_path, query])

# list of names for the headers
mash_df  = pd.read_csv(io.StringIO(mash_output.decode('utf-8')), sep='\t', header=None, names=['Reference', 'Query', 'distance', 'p-value', 'shared-hashes', 'ANI'])
number_hits = mash_df.shape[0]

#get the number of genomes wih mash distance < 0.2

if number_hits < 1:
    print(TerminalColors.RED +"""Error: No hits were found with the default settings
The phage likely represents a new species and genus 
However tax_my_phage is unable to classify it at this point as it can only classify at the 
          """)
    exit ()
else:
    print(TerminalColors.YELLOW +f"""
    Number of phage genomes detected with mash distance of < 0.2 is:{number_hits}""")

#sort dataframe by distance so they are at the top
mash_df = mash_df.sort_values(by='distance', ascending=True)
mash_df.to_csv('results/mash.txt', index=False)
minimum_value = mash_df['distance'].min()
maximum_value = mash_df.head(10)['distance'].max()

print(TerminalColors.BLUE +f"""The mash distances obtained for this query phage
is a minimum value of {minimum_value} and maximum value of {minimum_value} """)



#set the maximum number of hits to take forward. Max is 50 or the max number in the table if <50
filter_hits =""
if number_hits < 50:
    filter_hits = number_hits
else:
    filter_hits = 50

#copy  top 50 hits to a new dataframe
top_50 = mash_df.iloc[:filter_hits].copy()

ic (f"{mash_df.head(10)}")
ic (f"{top_50}")
#reindex 
top_50.reset_index(drop=True, inplace=True)

value_at_50th_position = top_50['distance'].iloc[filter_hits-1]
ic(f"{value_at_50th_position}")



#process each  row

for index, row in top_50.iterrows():
    dir_accesion = row.iloc[0]
    list_of_parts = dir_accesion.split('/')
    genus = list_of_parts[1]
    acc = list_of_parts[2]
    genus = genus.replace(".", "")
    acc = acc.replace(".fna", "")# swap fsa for fna in new ICTV.msh file 
    matches = taxa_df['Genbank'].str.contains(acc)
    matched_rows = taxa_df[matches].values.tolist()
    #result =matched_rows[0]
    id = row.iloc[2]
    ANI = (1-id)*100 
    subfamily = matched_rows[0][13]
    family  = matched_rows[0][12]
    genus  = matched_rows[0][14]
    top_50.loc[index, 'SubFamily'] = subfamily 
    top_50.loc[index, 'family'] = family
    top_50.loc[index, 'genus'] = genus
    top_50.loc[index, 'ANI'] = ANI


#what have you done below ? why was a dict being used ?

def get_genus_list (df):
   unique_genera = df.loc[:20, 'genus'].unique()
   ic(unique_genera)
   for genus in unique_genera:
       count_genera = (df.loc[:20, 'genus'] == genus).sum()
       ic(unique_genera)
       ic (count_genera)
   print (f"{count_genera}:{genus}")
   return unique_genera

#returns the unique genera names found in the mash hits - top_50 is not the best name!

unique_genera = get_genus_list(top_50)

#print for error checking
ic(unique_genera)

#number of genera
number_of_genera = len(unique_genera)

print(f"Found{number_of_genera} genera associated with this query genome")

#get all the keys for from a dictionary of accessions and genus names
keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]

#print the keys
ic(keys)


if len(unique_genera) == 1:
    print ("Only Found 1 genus so will proceed with getting all genomes associated with genus")
    keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]
    number_ok_keys =len(keys)
    print (f"Number of known species in the genus is {number_ok_keys} \n ")
    # create a command string for blastdbcmd
    get_genomes_cmd = f"blastdbcmd -db Bacteriophage_genomes.fasta  -entry {','.join(keys)} -out {known_taxa_path} "
    subprocess.run(get_genomes_cmd, shell=True, check=True)

elif len(unique_genera) >1:
    print("Found multiple genera that this query phage might be similar to so will proceed with processing them")
    list_of_genus_accessions =[]
    for i in unique_genera:
        keys = [k for k, v in accession_genus_dict.items() if v == i]
        number_of_keys = len(keys)
        #ic(keys)
        list_of_genus_accessions.extend(keys)
        print(f"Number of known species in the genus {i} is {number_of_keys}")
    ic(list_of_genus_accessions)
    ic(len(list_of_genus_accessions))
    get_genomes_cmd = f"blastdbcmd -db Bacteriophage_genomes.fasta  -entry {','.join(list_of_genus_accessions)} -out {known_taxa_path}"
    subprocess.run(get_genomes_cmd, shell=True, check=True)


def viridic_run(path_viridic_bash, results_path , viridic_in_path ):

    with open(known_taxa_path, 'r') as file1:
        with open(query, 'r') as file2:
            with open(viridic_in_path, 'w') as merged_file:
                merged_file.write(file1.read())
                merged_file.write(file2.read())
                #need to sort out the paths below still
                run_viridic = f" {path_viridic_bash} projdir={results_path} in={viridic_in_path} ncor={threads} sim_cols=Oranges  "
                subprocess.run(run_viridic, shell=True, check=True)

#get smallest mash distance
min_dist = top_50['distance'].min()

if  min_dist < 0.04: 
    print ("Phage is likely NOT a new species, will run further analysis now to to confirm this \n ")
    top_df = top_50[top_50['distance'] == min_dist]
    ic(top_df)


elif min_dist > 0.04 and min_dist  < 0.1: 
    print ("It is not clear if the phage is a new species or not. Will run further analysis now to confirm this...\n")
    top_df = top_50[top_50['distance'] < 0.1 ]
    ic(top_df)
    get_genus_list(top_df)

elif min_dist > 0.1  and min_dist  < 0.2:
    print ("Phage is a new species. Will run further analysis now ....\n")
    top_df = top_50[top_50['distance'] < 0.1 ]
    ic(top_df)


#######run viridic

viridic_run(path_viridic_bash, results_path , viridic_in_path )


def parse_virid():
    csv_file = os.path.join(results_path, "04_VIRIDIC_out","clusters.csv")
    df = pd.read_csv(csv_file, sep ="\t")
    ic(df.head)
    return df

#data fram from viridic created into a new dataframe df1
df1 = parse_virid()

taxa_df = pd.read_excel(VMR_path,sheet_name=0)

#Print the DataFrame
ic(taxa_df)
#renmae the columns again
taxa_df = taxa_df.rename(columns={'Virus GENBANK accession': 'Genbank'})
taxa_df['Genbank'].fillna('', inplace=True)

#merge the ICTV dataframe with the results of viridic
#fill in missing with Not Defined yet
merged_df = pd.merge(df1, taxa_df, left_on='genome', right_on='Genbank',how='left' ).fillna('Not Defined Yet')


#write dataframe to file
merged_df.to_csv(taxa_csv_output_path, sep='\t', index=False)

#create a copy of this dataframe for later use
copy_merged_df = merged_df.copy()


# Find the index of the row containing "taxmyPhage" in the "genome" column
index_to_remove = merged_df[merged_df['genome'] == 'taxmyPhage'].index
# Drop the row by index
merged_df.drop(index_to_remove, inplace=True)
# Reset the index to fill the gap created by removing the row
merged_df.reset_index(drop=True, inplace=True)


#Count the number genera
#excluding query
num_unique_viridic_genus_clusters = merged_df['genus_cluster'].nunique()
num_unique_ICTV_genera = merged_df['Genus'].nunique()


#including query
total_num_viridic_genus_clusters = copy_merged_df['genus_cluster'].nunique()
total_num_viridic_species_clusters = copy_merged_df['genus_cluster'].nunique()


print(f"""\n\nTotal number of VIRIDIC genus clusters in the input including QUERY sequence was:{total_num_viridic_genus_clusters}
Total number of VIRIDIC species clusters inluding QUERY sequence was {total_num_viridic_species_clusters} """)

print(f"""\n\nNumber of current ICTV defined genera was :{num_unique_ICTV_genera}
Number of VIRIDIC predicted genera (excluding query)was :{num_unique_viridic_genus_clusters} """)


if num_unique_ICTV_genera == num_unique_viridic_genus_clusters:
    print (f"""\n\nCurrent ICTV and VIRIDIC predictions are consistent for the data that was used to compare against""")




print(TerminalColors.BLUE +f"Number of unique VIRIDIC clusters at default cutoff of 70% is:{num_unique_viridic_genus_clusters}")
print(TerminalColors.BLUE +f"""Number of current ICTV genera associated with the reference genomes
 is {num_unique_ICTV_genera}""")

#unique_viridic_genus_clusters = merged_df['genus_cluster'].unique().tolist()
#num_unique_ICTV_genera = merged_df['Genus'].unique().tolist()

species_genus_dict = merged_df.set_index('species_cluster')['Species'].to_dict()



#get information on the query from the dataframe
#get species and genus cluster number
query_row = copy_merged_df[copy_merged_df['genome'] == 'taxmyPhage']
query_genus_cluster_number = query_row['genus_cluster'].values[0]
query_species_cluster_number = query_row['species_cluster'].values[0]

print(f"Cluster number of species is:{query_species_cluster_number} and cluster of genus is: {query_genus_cluster_number}")

print (f"Genus cluster number is {query_genus_cluster_number}  ")

#list of VIRIDIC genus and species numbers
list_ICTV_genus_clusters = merged_df['genus_cluster'].unique().tolist()
list_ICTV_species_clusters = merged_df['species_cluster'].unique().tolist()

ic (f"{list_ICTV_genus_clusters}")
ic (f"{list_ICTV_species_clusters}")

#create a dictionary linking genus_cluster to genus data
dict_genus_cluster_2_genus_name = merged_df.set_index('genus_cluster')['Genus'].to_dict()
dict_species_cluster_2_species_name = merged_df.set_index('species_cluster')['Species'].to_dict()

predicted_genus_name = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
print(f"Predicted genus is: {predicted_genus_name}")
#create a dict of species to species_cluster




#if number of ICTV genera and predicted VIRIDIC genera match:

if num_unique_ICTV_genera == num_unique_viridic_genus_clusters:
    print ("""Current ICTV taxonomy and VIRIDIC output appear to be consistent at the genus level""")

    #GENUS CHECK FIRST- Current genus and current species
    if query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number in list_ICTV_species_clusters:
        print(f"""Phage is within a current genus and same as a current species 
         ....working out which one now .....""")
        predicted_genus = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
        predicted_species = dict_species_cluster_2_species_name[query_species_cluster_number]
        print(f"""QUERY is in the genus:{predicted_genus} and is species: {predicted_species}""")
        #identify the row in the pandas data frame that is the same species
        matching_species_row = merged_df[merged_df['Species'] == predicted_species]
        ic (f"{matching_species_row}")
        list_of_S_data = matching_species_row[0:].values.flatten().tolist()
        ic(f"{list_of_S_data[14:20]}")
        ic(f"{list_of_S_data}")
        print(TerminalColors.YELLOW + f"""Query sequence is: 
                Class:{list_of_S_data[14]}
                Family: {list_of_S_data[15]}
                Subfamily:{list_of_S_data[16]}
                Genus:{list_of_S_data[17]}
                Species:{list_of_S_data[19]}
                 """)
        with open(summary_output_path,'a') as file:
            file.write(f"""
Query sequence can be classified within a current genus and species, it is in:\n
Class:{list_of_S_data[14]}\tFamily: {list_of_S_data[15]}\tSubfamily:{list_of_S_data[16]}\tGenus:{list_of_S_data[17]}Species:{list_of_S_data[19]}
\n The data from the initial mash searching is below as tsv format \n
Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing determine if it represents a new 
species or geneus. It does not tell you if it is similar to other phages that have yet to be classified
You can do this by comparison with INPHARED database if you wish""")
        mash_df.to_csv(summary_output_path, mode='a', header=True, index=False,sep='\t')

        #WRITE CODE FOR GIVING INFO ON SPECIES

    #SAME GENUS but different species
    elif query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number not in list_ICTV_species_clusters:
        print(f"""Phage is within a current genus, BUT is representative of a new species 
                 ....working out which one now .....""")

        matching_genus_rows = merged_df[merged_df['genus_cluster'] == query_genus_cluster_number]
        dict_exemplar_genus = matching_genus_rows.iloc[0].to_dict()
        genus_value =dict_exemplar_genus['Genus']
        ic(f"f{matching_genus_rows}")
        ic(f"{genus_value}")

        print(TerminalColors.YELLOW + f"""Query sequence is: 
        Class:{dict_exemplar_genus['Class']}
        Family: {dict_exemplar_genus['Family']} 
        Subfamily:{dict_exemplar_genus['Subfamily']}
        Genus:{dict_exemplar_genus['Genus']}
        Species:{dict_exemplar_genus['Genus']} new_name
         """)

        with open(summary_output_path, 'a') as file:
            file.write(f"""
        Query sequence can be classified within a current genus and represents a new species, it is in:\n
        Class:{dict_exemplar_genus['Class']}\tFamily: {dict_exemplar_genus['Family']}\tSubfamily:{dict_exemplar_genus['Subfamily']}\tGenus:{dict_exemplar_genus['Genus']}\tnew_specices_name
        \n The data from the initial mash searching is below as tsv format \n
        Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing determine if it represents a new 
        species or geneus. It does not tell you if it is similar to other phages that have yet to be classified
        You can do this by comparison with INPHARED database if you wish""")
        mash_df.to_csv(summary_output_path, mode='a', header=True, index=False,sep='\t')



    elif query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number not in list_ICTV_species_clusters:
        print(f"""Query does not fall within  a  current genus or species as defined by ICTV
        Therefore the query sequence is likely the first representative of both a new species and new genus
        Data produced by taxmyPHAGE will help you write a Taxonomy proposal so it can be offically classified
        WARNING taxmyPHAGE does not compare against all other known phages, only those that have been classified
        """)

        with open(summary_output_path, 'a') as file:
            file.write(f"""
        Query sequence can not be classified within a current genus or species, it is in:\n
        Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing determine if it represents a new 
        species or geneus. It does not tell you if it is similar to other phages that have yet to be classified
        You can do this by comparison with INPHARED database if you wish""")
        mash_df.to_csv(summary_output_path, mode='a', header=True, index=False,sep='\t')



######if number of VIRIDIC genera is greater than ICTV genera
elif num_unique_ICTV_genera < num_unique_viridic_genus_clusters:
    print (TerminalColors.RED+"""The number of expected genera based on current ICTV classification is less than the predicted 
number of genus clusters as predicted by VIRIDIC.
This does not mean the current ICTV classification is wrong (it might be) or that VIRIDIC
is wrong. It could be an edge case that automated process cannot distinguish. It will require
 more manual curation to look at the output files """)
    if query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number in list_ICTV_species_clusters:
        print(TerminalColors.BLUE+"""Phage is within a current genus and same as a current species 
         ....working out which one now .....""")
        predicted_genus = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
        predicted_species = dict_species_cluster_2_species_name[query_species_cluster_number]
        print(f"""QUERY is in the genus:{predicted_genus} and is species: {predicted_species}""")
        #identify the row in the pandas data frame that is the same species
        matching_species_row = merged_df[merged_df['Species'] == predicted_species]
        ic (f"{matching_species_row}")
        list_of_S_data = matching_species_row[0:].values.flatten().tolist()
        ic(f"{list_of_S_data[14:20]}")
        ic(f"{list_of_S_data}")
        print(TerminalColors.YELLOW + f"""Query sequence is in the;
                Class:{list_of_S_data[14]}
                Family: {list_of_S_data[15]}
                Subfamily:{list_of_S_data[16]}
                Genus:{list_of_S_data[17]}
                Species:{list_of_S_data[19]}
                 """)

        with open(summary_output_path, 'a') as file:
            file.write(f"""\n
The number of expected genera based on current ICTV classification is less than the predicted 
number of genus clusters as predicted by VIRIDIC.
This does not mean the current ICTV classification is wrong (it might be) or that VIRIDIC
is wrong. It could be an edge case that automated process cannot distinguish. It will require
more manual curation to look at the output files\n 
Query sequence can be classified within a current genus and species, it is in:\n
Class:{dict_exemplar_genus['Class']}\tFamily: {dict_exemplar_genus['Family']}\tSubfamily:{dict_exemplar_genus['Genus']}\tGenus:{list_of_S_data[17]}\tnew_specices_name
\n The data from the initial mash searching is below as tsv format \n
Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing you determine if it represents a new 
species or genus. It does not tell you if it is similar to other phages that have yet to be classified
You can do this by comparison with INPHARED database if you wish\n\n""")
        mash_df.to_csv(summary_output_path, mode='a', header=True, index=False,sep='\t')



    elif query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number not in list_ICTV_species_clusters:
        print(TerminalColors.BLUE+"""Phage is within a current genus, BUT is representative of a new species 
                 ....working out which one now .....""")

        matching_genus_rows = merged_df[merged_df['genus_cluster'] == query_genus_cluster_number]
        dict_exemplar_genus = matching_genus_rows.iloc[0].to_dict()
        genus_value =dict_exemplar_genus['Genus']
        ic(f"f{matching_genus_rows}")
        ic(f"{genus_value}")

        print(TerminalColors.YELLOW + f"""Query sequence is in the;
        Class:{dict_exemplar_genus['Class']}
        Family: {dict_exemplar_genus['Family']} 
        Subfamily:{dict_exemplar_genus['Subfamily']}
        Genus:{dict_exemplar_genus['Genus']}
        Species:new_species_name
         """)

        with open(summary_output_path, 'a') as file:
            file.write(f"""\n
Query sequence can be classified within a current genus and represents a new species, it is in:\n
Class:{dict_exemplar_genus['Class']}\tFamily: {dict_exemplar_genus['Family']}\tSubfamily:{dict_exemplar_genus['Subfamily']}\tGenus:{dict_exemplar_genus['Genus']}\tSpecies:name_your_species
\n 
Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing you to determine if it represents a new 
species or genuss. It does not tell you if it is similar to other phages that have yet to be classified
You can do this by comparison with INPHARED database if you wish\n
The data from the initial mash searching is below as tsv format \n\n""")
        mash_df.to_csv(summary_output_path, mode='a', header=True, index=False,sep='\t')


#clean up files that are left over
dirs_to_remove = ['01_BlastDB', '02_BlastN_out', '03_calculations_out']

for dir in dirs_to_remove:
    to_remove = os.path.join(results_path, dir)
    ic(f"{to_remove}")
    if os.path.isdir(to_remove):
        shutil.rmtree(to_remove)




#copy and rename heatmap files

heatmap_in = os.path.join(results_path,'04_VIRIDIC_out','Heatmap.PDF')
heatmap_out =os.path.join(results_path, 'viridic_heatmap.pdf')

shutil.copy(heatmap_in, heatmap_out)

# Path to the input PDF file
input_pdf_path = 'viridic_heatmap.pdf'

# Path to the output SVG file
output_svg_path = os.path.join(results_path, 'output.svg')

convert_pdf_cmd = f" pdf2svg {heatmap_out} {output_svg_path}"
# execute the command using subprocess
subprocess.run(convert_pdf_cmd, shell=True, check=True)

#keep 04 folder if verbose

if not verbose:
    to_remove = os.path.join(results_path,"04_VIRIDIC_out")
    if os.path.isdir(to_remove):
        shutil.rmtree(to_remove)





# #print (merged_df.head(10))
#
# #create dictionary of key value pairs
# dict_genome = merged_df.set_index('genome')[['species_cluster', 'genus_cluster','Species','Genus','Family']].to_dict('index')
#
# #get a list of unique  genera clusters
# unique_genera = merged_df['genus_cluster'].unique().tolist()
#
# #get a list of unique genera names
# unique_genera_names = merged_df['Genus'].unique().tolist()
#
# #get number of genus clusters
# num_genera=len(unique_genera)
#
# #get number of genus names
# num_genus_names=len(unique_genera_names)
#
# #get a list of unquie species clusters
# unique_species_cluster = merged_df['species_cluster'].unique().tolist()
#
# #get a list of unique species  names
# unique_species_names  = merged_df['Species'].unique().tolist()
#
# #get number of species clusters
# num_species=len(unique_species_cluster)
#
# #get number of species names
#
# ic(unique_species_names)
#
# ic (f"Number of species is:{num_species}")
#
# ic(num_genus_names)
#
# ic (f"Number of genera is:{num_genera}")
#
# if num_genera+1 > num_genus_names:
#      print("There could be a problem with classification of this phage. Further manual checking \
# of the results might be required")
#
# #get a list of lists for the query genome
# query_row_list  = merged_df.loc[merged_df['genome'] == 'taxmyPhage'].values.tolist()
#
# #extract the "species" cluster number from the list of lists. It is the 2nd element of the 1st list
# query_cluster_number = query_row_list[0][1]
#
# count_query_cluster_in_df  = merged_df['species_cluster'].value_counts()[query_cluster_number]
#
# if count_query_cluster_in_df == 1:
#     print (f"The phage represents a novel species")
# elif count_query_cluster_in_df != 1 :
#     print (f"The phage is the same as an existing species")
#
# ic(count_query_cluster_in_df)
#
# ic(query_cluster_number)
#

# shutil.copy('04_VIRIDIC_out/Heatmap.PDF', './viridic_heatmap.pdf')
#
# # Path to the input PDF file
# input_pdf_path = 'viridic_heatmap.pdf'
#
# # Path to the output SVG file
# output_svg_path = 'output.svg'
#
# convert_pdf_cmd = f" pdf2svg {input_pdf_path} {output_svg_path}"
# # execute the command using subprocess
# subprocess.run(convert_pdf_cmd, shell=True, check=True)
#
# #still do to
# #cleanup files that are left over
#
# #merged_df = merged_df.loc[merged_df['Genus'] != 'Not Defined Yet']
#
# grouped = merged_df.groupby('Genus')['genus_cluster'].nunique()
# same_cluster_numbers = len(grouped.unique()) == 1
#
# if same_cluster_numbers:
#     print("All genera have the same genus_cluster numbers.")
# else:
#     print("Not all genera have the same genus_cluster numbers.")
#
# #group genus and genus_cluster
# grouped = merged_df.groupby(['Genus']).agg({'genus_cluster': 'nunique'})
#
# # Filter to show only the rows where the number of unique genus_cluster values is greater than 1
# nonuniform_clusters = grouped[grouped['genus_cluster'] >= 1]
#
# # Print the resulting rows
# if len(nonuniform_clusters) > 0:
#     print("The following genera have non-uniform genus_cluster numbers:")
#     print(nonuniform_clusters)
#     print ("This is inconsistent with automated taxonomy as determined using VIRIDC for cut offs \
# there are a number of reasons for this. It could be as a result of incorrect allocation of \
# species to a genera, the genus had not recently been updated with current cutoffs or real \
# differences. To understand this further you will need to read the taxonomy proposal on the ICTV website ")
# else:
#     print("All genera have the same genus_cluster number")
#
# genus_cluster_matching_query = merged_df.loc[merged_df['genus_cluster'] == query_cluster_number ]
#
# ic(genus_cluster_matching_query)
#
# taxmyPhage_clusters = merged_df.loc[merged_df['genome'] == 'taxmyPhage', ['species_cluster', 'genus_cluster']].values.tolist()[0]
#
# ic (taxmyPhage_clusters)
