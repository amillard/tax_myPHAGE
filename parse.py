import pandas as pd 
import shutil
import subprocess


from icecream import ic

def parse_virid():
    csv_file ="04_VIRIDIC_out/clusters.csv"
    df = pd.read_csv(csv_file, sep ="\t")
    print(df.head)
    return df
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
     print("There could be a problem with classificaiton of this phage. Further manual checking \
of the results might be required")

#get a list of lists for the query genome 
query_row_list  = merged_df.loc[merged_df['genome'] == 'UP30'].values.tolist()

#extract the "species" cluster number from the list of lists. It is the 2nd element of the 1st list 
query_cluster_number_is = query_row_list[0][1]

count_query_cluster_in_df  = merged_df['species_cluster'].value_counts()[query_cluster_number_is]

if count_query_cluster_in_df == 1:
    print (f"The phage represents a novel species")

ic(count_query_cluster_in_df)

ic(query_cluster_number_is)

dirs_to_remove = ['01_BlastDB', '02_BlastN_out', '03_calculations_out']

#for directory in dirs_to_remove:
#    shutil.rmtree(directory)

shutil.copy('04_VIRIDIC_out/Heatmap.PDF', './viridic_heatmap.pdf')

import os
from pdf2image import convert_from_path
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM

# Path to the input PDF file
input_pdf_path = 'viridic_heatmap.pdf'

# Path to the output SVG file
output_svg_path = 'output.svg'

convert_pdf_cmd = f" pdf2svg {input_pdf_path} {output_svg_path}"
# execute the command using subprocess
subprocess.run(convert_pdf_cmd, shell=True, check=True)


