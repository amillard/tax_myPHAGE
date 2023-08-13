import pandas as pd
import shutil
import sys
import os 
from icecream import ic

#requires GenomesDB from INPHARED to build the database

# this is the location of where the script and the databases are (instead of current_directory which is the users current directory)
HOME = os.path.dirname(__file__)


#check the files are present
VMR_path = os.path.join(HOME, 'VMR.xlsx')


#get a datarame
taxa_df = pd.read_excel(VMR_path,sheet_name=0)

#filter on host bacteria and entries that are present in Genbank
phage_df = taxa_df[(taxa_df['Host source'] == 'bacteria') & (taxa_df['Genome coverage'] != 'No entry in Genbank' )]

#check the size
print (phage_df.shape)

#filter out pesky segments
phage_df = phage_df[~phage_df['Virus GENBANK accession'].str.contains(';')]
print (phage_df.shape )

#drop empty genera and get a list of empty genera
cp_df = phage_df[phage_df['Genus'].isna()].copy
ic(cp_df)
phage_df = phage_df.dropna(subset=['Genus'])

acc = phage_df['Virus GENBANK accession'].tolist()

#specify the GENOMESDB Folder
genomesDB = 'GenomesDB'

#create a column with the extra data or path to the file
phage_df['checkpath']=  phage_df['Genus']+ '/' + phage_df['Virus GENBANK accession'] + '.fna'


list = phage_df['checkpath'].tolist()


#

not_found = 0
for i in list:
    #print(f"Searhcing for:Genera \t {i}")
    filepath = os.path.join(HOME,'Genera',i)

    if os.path.isfile(filepath):
       pass
       # print (f"found {filepath}")
    else:
        print  (f"Warning {filepath} not FOUND !!")
        not_found +=1

#check the msh file is present in each subfolder

for i in list:
    msh = i+'.msh'
    filepath = os.path.join(HOME, 'Genera', msh)
    if os.path.isfile(filepath):
       pass
       # print (f"found {filepath}")
    else:
        print  (f"Warning {filepath} not FOUND !!")
        not_found +=1


