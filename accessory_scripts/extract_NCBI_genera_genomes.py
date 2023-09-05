#Link to INPHARED data to get all genomes that NCBI classify within a genus
#Uses the data from INPHARED Conversion currently, so have to add the  column headers to the df

import pandas as pd
import subprocess
import os
from icecream import ic
from argparse import ArgumentParser

usage = "%prog [options] file (or - for stdin)"
description= """Extract all genomes classified by GenBank at the taxa level to a single file
Note this is not the same as ICTV defined genomes  """
parser = ArgumentParser(usage, description=description)
parser.add_argument("-v", "--verbose", action="store_true", default = 0)
parser.add_argument('-g', '--genus', dest='genus' , help="you will need to know the genus name in advance", type=str)
parser.add_argument("-df", dest='datafile', type=str, help="Specify full path to a INPHARED monthly data file "
                                                 "eg ./25August2023_data.tsv ")


args, nargs = parser.parse_known_args()
verbose = args.verbose

#turn on ICECREAM reporting
if not verbose: ic.disable()


data_file = args.datafile
genus = args.genus
print (f"{data_file} {genus}")


df = pd.read_csv(data_file, sep='\t', header=None, names=['Accession','Date Updated','Genetic Material','Phage Description','Length','GC content','Realm','Kingdon','Phylum','Class','Order','Sub-family','Family','Genus','Baltimore Group','JumboPhage','Coding Capicity','Host','Host2'])

genus_df = df[df["Genus"] == genus]
genus_df = genus_df[~genus_df['Accession'].str.contains("NC_")]

#filter genomes out that are very small compared to median
median_genome_length = genus_df['Length'].median()
std_dev_genome_length = genus_df['Length'].std()

min_genome_threshold = median_genome_length - (median_genome_length/10)

genus_df = genus_df[genus_df["Length"] >= min_genome_threshold ]

print(f"{genus_df}, {median_genome_length},{std_dev_genome_length},{min_genome_threshold}")

list_of_genus_acc = genus_df["Accession"].tolist()

print (f"{list_of_genus_acc}")


input ="1Aug2023_genomes.fa"

makeblastdb_command = f"makeblastdb -in {input} -parse_seqids -dbtype nucl"

ndb_file = input+".ndb"

if os.path.isfile(ndb_file):
    print ("BlastDB already exists will skip this step")
else:
    try:
        subprocess.run(makeblastdb_command, shell=True, check=True)
        print("makeblastdb command executed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing makeblastdb: {e}")


def extract_genomes(blastdb,list_acc,genus):
    """

    :param blastdb:
    :param list_acc:
    :param genus:
    :param output:
    :return:
    """
    output = genus+"_all_genomes.fna"
    get_genomes_cmd = f"blastdbcmd -db {blastdb}  -entry {','.join(list_acc)} -out {output} "
    res = subprocess.getoutput(get_genomes_cmd)

extract_genomes(input,list_of_genus_acc,genus)

