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



# Create the argument parser
usage = "%prog [options] file (or - for stdin)"
description='Read input.fasta file'
parser = ArgumentParser(usage, description=description)
parser.add_argument("-v", "--verbose", action="store_true", default = 0)
parser.add_argument('-i', '--input', dest='in_fasta', type=str, help='Path to input.fasta file')
args, nargs = parser.parse_known_args()
verbose = args.verbose
fasta_file = args.in_fasta

if not verbose: ic.disable()


query ="query.fasta"

#write the input file to a new file with a header called "taxmyPhage" 
handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file)
entries = SimpleFastaParser(handle)

with open(args.in_fasta, "r") as input_file:
    # Open output FASTA file
    with open(query, "w") as output_file:
        # Loop over input sequences
        for name, seq in entries:
            print(f">taxmyPhage\n{seq}", file=output_file)


current_directory = os.getcwd()

blastdb_path = os.path.join(current_directory, 'Bacteriophage_genomes.fasta')


#check the VMR is present 
VMR_path = "VMR.xls"

if os.path.exists(VMR_path):
    print (f"Found {VMR_path} as expected")
else:
    print(f'File {VMR_path} does not exist')

#check if blastDB is present
if os.path.exists(VMR_path):
    print (f"Found {blastdb_path} as expected")
else:
    print(f"File {blastdb_path} does not exist will create database now  ")
    print("Still do do this bit of installing the databases automatically")
    



#Read the viral master species record into a DataFrame
taxa_df = pd.read_excel('VMR.xls',sheet_name=0)

# Print the DataFrame
print(taxa_df)
taxa_df = taxa_df.rename(columns={'Virus GENBANK accession': 'Genbank'})
taxa_df['Genbank'].fillna('', inplace=True)
# Get the column headings as a list

headings = list(taxa_df.columns)
#crate a dictionary of Accessions linking to Genus 
accession_genus_dict = taxa_df.set_index('Genbank')['Genus'].to_dict()

# Print the column headings one per line
for heading in headings:
    ic(heading)


#run mash to get top hit and read into a pandas dataframe

mash_output = subprocess.check_output(['mash', 'dist', '-d' , '0.25',  '-p', '8', 'ICTV.msh', query])

# list of names for the headers
mash_df  = pd.read_csv(io.StringIO(mash_output.decode('utf-8')), sep='\t', header=None, names=['Reference', 'Query', 'distance', 'p-value', 'shared-hashes', 'ANI'])

#sort dataframe to get  top hits first 
mash_df  = mash_df.sort_values(by='distance', ascending=True)

ic(mash_df)
ic(mash_df.tail())
#copy  top 50 hits to a new dataframe
top_50 = mash_df.iloc[:150].copy()
#reindex 
top_50.reset_index(drop=True, inplace=True)
ic(top_50)

#process each  row 

for index, row in top_50.iterrows():
    dir_accesion = row.iloc[0]
    list_of_parts = dir_accesion.split('/')
    genus = list_of_parts[1]
    acc = list_of_parts[2]
    genus = genus.replace(".", "")
    acc = acc.replace(".fsa", "")
#   print (f"{acc}")
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
      # ic(unique_genera)
      # ic (count_genera)
   print (f"{count_genera}:{genus}")
   return unique_genera

#returns the unique genera names found in the mash hits
unique_genera = get_genus_list(top_50)

ic(unique_genera)

number_of_genera = len(unique_genera)

print (f"Found:{number_of_genera} genera associated with this query genome")

keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]
#ic(keys)
####

if len(unique_genera) == 1:
    print ("Only Found 1 genus so will proceed with getting all genomes associated with genus")
    keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]
    number_ok_keys =len(keys)
    print (f"Number of known genomes in the genus is {number_ok_keys}")
    # create a command string for blastdbcmd
    get_genomes_cmd = f"blastdbcmd -db Bacteriophage_genomes.fasta  -entry {','.join(keys)} -out known_taxa.fa"
    subprocess.run(get_genomes_cmd, shell=True, check=True)

elif len(unique_genera) >1:
    print("Found multiple genera that this query phage might be similar to so will proceed with processing them")
    list_of_genus_accessions =[]
    for i in unique_genera:
        keys = [k for k, v in accession_genus_dict.items() if v == i]
        number_of_keys = len(keys)
        #ic(keys)
        list_of_genus_accessions.extend(keys)
        print(f"Number of known genomes in the genus is {number_of_keys}")
    #ic(list_of_genus_accessions)
    ic (len(list_of_genus_accessions))
    get_genomes_cmd = f"blastdbcmd -db Bacteriophage_genomes.fasta  -entry {','.join(list_of_genus_accessions)} -out known_taxa.fa"
    subprocess.run(get_genomes_cmd, shell=True, check=True)


def viridic_run():
    with open('known_taxa.fa', 'r') as file1:
        with open(query, 'r') as file2:
            with open('viridic_in.fa', 'w') as merged_file:
                merged_file.write(file1.read())
                merged_file.write(file2.read())
                run_viridic = f" ~/viridic_v1.0_r3.6/viridic.bash projdir=/home/andy/taxa_my_phage/ in=/home/andy/taxa_my_phage/viridic_in.fa"
                #subprocess.run(run_viridic, shell=True, check=True)

#get smallest mash distance
min_dist = top_50['distance'].min()


if  min_dist < 0.04: 
    print ("Phage is likely NOT a new species, will run futher analysis now to to confirm this")
    top_df = top_50[top_50['distance'] == min_dist]
    print(top_df) 


elif min_dist > 0.04 and min_dist  < 0.1: 
    print ("It is not clear if the phage is a new species or not. Will run further analysis now to confirm this")
    top_df = top_50[top_50['distance'] < 0.1 ]
    print(top_df)
    get_genus_list(top_df,accession_genus_dict)

elif min_dist > 0.1  and min_dist  < 0.2:
    print ("Phage is a new species. Will run further analysis now   ")
    top_df = top_50[top_50['distance'] < 0.1 ]
    print(top_df)


#run viridic based on above
viridic_run()



def parse_virid():
    csv_file ="04_VIRIDIC_out/clusters.csv"
    df = pd.read_csv(csv_file, sep ="\t")
    print(df.head)
    return df

#data fram from viridic
df1 = parse_virid()

taxa_df = pd.read_excel('VMR.xls',sheet_name=0)

# Print the DataFrame
#print(taxa_df)
taxa_df = taxa_df.rename(columns={'Virus GENBANK accession': 'Genbank'})
taxa_df['Genbank'].fillna('', inplace=True)

#merge the ICTV dataframe with the results of viridic
merged_df = pd.merge(df1, taxa_df, left_on='genome', right_on='Genbank',how='left' ).fillna('Not Defined Yet')

#write dataframe to file
merged_df.to_csv('Taxa_Output.csv', sep='\t', index=False)

#create dictionary of key value pairs
dict_genome = merged_df.set_index('genome')[['species_cluster', 'genus_cluster','Species','Genus','Family']].to_dict('index')

#get a list of unique  genera clusters
unique_genera = merged_df['genus_cluster'].unique().tolist()

#get a list of unique genera names
unique_genera_names = merged_df['Genus'].unique().tolist()

#get number of genus clusters
num_genera=len(unique_genera)

#get number of genus names
num_genus_names=len(unique_genera_names)

#get a list of unquie species clusters
unique_species_cluster = merged_df['species_cluster'].unique().tolist()

#get a list of unique species  names
unique_species_names  = merged_df['Species'].unique().tolist()

#get number of species clusters
num_species=len(unique_species_cluster)

#get number of species names

ic(unique_species_names)

ic (f"Number of species is:{num_species}")

ic(num_genus_names)

ic (f"Number of genera is:{num_genera}")

if num_genera+1 > num_genus_names:
     print("There could be a problem with classification of this phage. Further manual checking \
of the results might be required")

#get a list of lists for the query genome
query_row_list  = merged_df.loc[merged_df['genome'] == 'taxmyPhage'].values.tolist()

#extract the "species" cluster number from the list of lists. It is the 2nd element of the 1st list
query_cluster_number = query_row_list[0][1]

count_query_cluster_in_df  = merged_df['species_cluster'].value_counts()[query_cluster_number]

if count_query_cluster_in_df == 1:
    print (f"The phage represents a novel species")
elif count_query_cluster_in_df != 1 :
    print (f"The phage is the same as an existing species")

ic(count_query_cluster_in_df)

ic(query_cluster_number)

dirs_to_remove = ['01_BlastDB', '02_BlastN_out', '03_calculations_out']

#for directory in dirs_to_remove:
for dir in dirs_to_remove:
    if os.path.isdir(dir):
        shutil.rmtree(dir)

shutil.copy('04_VIRIDIC_out/Heatmap.PDF', './viridic_heatmap.pdf')

# Path to the input PDF file
input_pdf_path = 'viridic_heatmap.pdf'

# Path to the output SVG file
output_svg_path = 'output.svg'

convert_pdf_cmd = f" pdf2svg {input_pdf_path} {output_svg_path}"
# execute the command using subprocess
subprocess.run(convert_pdf_cmd, shell=True, check=True)

#still do to 
#cleanup files that are left over 

#merged_df = merged_df.loc[merged_df['Genus'] != 'Not Defined Yet']

grouped = merged_df.groupby('Genus')['genus_cluster'].nunique()
same_cluster_numbers = len(grouped.unique()) == 1

if same_cluster_numbers:
    print("All genera have the same genus_cluster numbers.")
else:
    print("Not all genera have the same genus_cluster numbers.")

#group genus and genus_cluster
grouped = merged_df.groupby(['Genus']).agg({'genus_cluster': 'nunique'})

# Filter to show only the rows where the number of unique genus_cluster values is greater than 1
nonuniform_clusters = grouped[grouped['genus_cluster'] > 1]

# Print the resulting rows
if len(nonuniform_clusters) > 0:
    print("The following genera have non-uniform genus_cluster numbers:")
    print(nonuniform_clusters)
    print ("This is inconsistent with automated taxonomy as determined using VIRIDC for cut offs \
there are a number of reasons for this. It could be as a result of incorrect allocation of \
species to a genera, the genus had not recently been updated with current cutoffs or real \
differences. To understand this further you will need to read the taxonomy proposal on the ICTV website ")
else:
    print("All genera have the same genus_cluster number.")

genus_cluster_matching_query = merged_df.loc[merged_df['genus_cluster'] == query_cluster_number ]

ic(genus_cluster_matching_query)

#taxmyPhage_clusters = merged_df.loc[merged_df['genome'] == 'taxmyPhage', ['species_cluster', 'genus_cluster']].values.tolist()[0]

#ic (taxmyPhage_clusters)
