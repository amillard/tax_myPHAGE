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
from itertools import combinations_with_replacement
import networkx as nx
from subprocess import getoutput
from tqdm import tqdm
import time
from datetime import timedelta
import random
import matplotlib.pyplot as plt
import numpy as np


def print_error(txt):
    print(f"\033[31m{txt}\033[0m")
def print_warn(txt):
    print(f"\033[94m{txt}\033[0m")
def print_ok(txt):
    print(f"\033[34m{txt}\033[0m")
def print_res(txt):
    print(f"\033[33m{txt}\033[0m")

class PoorMansViridic:
    def __init__(self, file, genus_threshold=70, species_threshold=95, nthreads=1, verbose = True):
        self.verbose = verbose
        self.file = file
        self.result_dir = os.path.dirname(self.file)
        self.nthreads = nthreads
        self.genus_threshold = genus_threshold
        self.species_threshold = species_threshold

    def run(self):
        self.makeblastdb()
        self.blastn()
        self.parse_blastn_file()
        self.calculate_distances()
        self.cluster_all()
        return self.dfT, self.pmv_outfile

    def cluster_all(self):
        dfTg = self.sim2cluster(self.genus_threshold, 'genus')#.rename({'cluster':'genus_cluster'}, axis=1)
        dfTs = self.sim2cluster(self.species_threshold, 'species')#.rename({'cluster':'species_cluster'}, axis=1)
        dfT = pd.merge(dfTg, dfTs, on='genome').sort_values('species_cluster genus_cluster'.split())
        dfT.reset_index(drop=True, inplace=True)
        self.pmv_outfile = os.path.join(self.result_dir, os.path.basename(self.file) + '.genus_species_clusters.tsv')
        dfT.to_csv(self.pmv_outfile, index=False, sep='\t')
        self.dfT = dfT

    def sim2cluster(self, th, tax_level):
        ic(f"Generating graph for finding", tax_level, "clusters")
        M = self.dfM
        G = nx.from_pandas_edgelist(M[(M.sim >= th) & (M.A != M.B)], source='A', target='B')
        singletons = list(set(M.A.unique().tolist()).difference(G.nodes()))
        G.add_nodes_from(singletons)

        graphs = [G.subgraph(x) for x in nx.connected_components(G)]
        L = []
        for n, g in enumerate(graphs):
            L.extend([(node, n+1) for node in g.nodes()])

        return pd.DataFrame(L, columns=f'genome {tax_level}_cluster'.split())

    def makeblastdb(self):
        outfile = self.file + '.nin'
        if not os.path.exists(outfile):
            cmd = f'makeblastdb -in {self.file}  -dbtype nucl'
            ic("Creating blastn database:", cmd)
            res = getoutput(cmd)
            print(res)
        
    def blastn(self):
        outfile = os.path.join(self.result_dir, os.path.basename(self.file) + '.blastn_vs_self.tab')
        if not os.path.exists(outfile):
            cmd = f'blastn -evalue 1 -max_target_seqs 10000 -num_threads {self.nthreads} -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {self.file} -db {self.file} -outfmt "6 qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore" -out {outfile}'
            ic("Blasting against itself:", cmd)
            print(cmd)
            getoutput(cmd)

        self.blastn_result_file = outfile

    def parse_blastn_file(self):
        ic("Reading", self.blastn_result_file)
        df = pd.read_table(self.blastn_result_file, names='qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore'.split())

        self.size_dict  = dict(df['qseqid qlen'.split()].values) | dict(df['sseqid slen'.split()].values)

        # this part is slow ... filling the genome with all the matches
        M = {}
        for qseqid, sseqid, pident, length, qlen, slen, mismatch, nident, gapopen, qstart, qend, sstart, send, qseq, sseq, evalue, bitscore in tqdm(df.values):
            key = (qseqid, sseqid)
            if qstart > qend: qstart, qend = qend, qstart
            M.setdefault(key, [0]*qlen)
            v = M[key]
            gaps = 0
            for n, nt in enumerate(qseq):
                if nt == '-':
                    gaps += 1
                    continue
                if v[qstart-1+n-gaps] == 1: continue
                v[qstart-1+n-gaps] = int(sseq[n] == qseq[n])

        self.M = M

    def calculate_distances(self):
        M = self.M
        size_dict = self.size_dict
        genomes = list(set([x[0] for x in M.keys()] + [x[0] for x in M.keys()]))
        # now calculating the distances 
        L = []
        for gA, gB in combinations_with_replacement(genomes, 2):
            key = (gB, gA)
            v = M[key]
            idAB = sum(M[(gA, gB)])
            idBA = sum(M[(gB, gA)])
            lA = size_dict[gA]
            lB = size_dict[gB]
            simAB = ((idAB+idBA)*100)/(lA+lB)
            distAB = 100-simAB

            aligned_fraction_genome_1 = idAB / lA
            aligned_fraction_genome_2 = idBA / lB
            genome_length_ratio = lA / lB if lA < lB else lB / lA
            L.append([gA, gB, distAB, aligned_fraction_genome_1, genome_length_ratio, aligned_fraction_genome_2])

        dfM = pd.DataFrame(L, columns='A B distAB afg1 glr afg2'.split())

        dfM['sim'] = 100 - dfM.distAB
        self.dfM = dfM

    def save_similarities(self, outfile='similarities.tsv'):
        df = self.dfM['A B sim'.split()]
        df = df[df.A != df.B]
        df.sort_values('sim', ascending=False)
        df.index.name = ''
        df.to_csv(outfile, index=False, sep='\t')
        
    
def heatmap(dfM, outfile, cmap='cividis'):

    svg_out = outfile+".svg"
    pdf_out = outfile+".pdf"
    jpg_out = outfile+".jpg"
    ax = plt.gca()
    dfM.update(dfM.loc[dfM.A > dfM.B].rename({'A': 'B', 'B': 'A'}, axis=1))
    dfM = dfM.round(2)
    df = dfM.pivot(index='A', columns='B', values='sim').fillna(0)
    df = df.rename({'taxmyPhage':'query'}, axis=1).rename({'taxmyPhage':'query'}, axis=0)
    im = plt.imshow(df.values, cmap=cmap)
    ax.set_xticks(np.arange(df.shape[1]), labels=df.columns.tolist())
    ax.set_yticks(np.arange(df.shape[0]), labels=df.index.tolist())

    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(df.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(df.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            text = ax.text(j, i, df.iloc[i, j], ha="center", va="center", color="w")

    plt.savefig(svg_out)
    plt.savefig(pdf_out)
    plt.savefig(jpg_out)
    return im

        
def is_program_installed_unix(program_name):
    try:
        subprocess.check_output(f'which {program_name}', shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def check_programs():
    # check programs are installed
    program_name = "blastdbcmd"
    if is_program_installed_unix(program_name):
        ic(f"{program_name} is installed will proceed ")
    else:
        print_error(f"{program_name} is not installed.")
        exit()

    program_name = "mash"
    if is_program_installed_unix(program_name):
        ic(f"{program_name} is installed will proceed ")
    else:
        print_error(f"{program_name} is not installed.")
        exit()

def check_blastDB(blastdb_path):
    #check if blastDB is present
    if os.path.exists(blastdb_path):
        print_ok(f"Found {blastdb_path} as expected")
    else:
        print_error(f"File {blastdb_path} does not exist will create database now  ")
        print_error("Will download the database now and create database")
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

    
if __name__ == '__main__': 
    usage = "%prog [options] file (or - for stdin)"
    description= """Takes a phage genome as as fasta file and compares against all phage genomes that are currently classified 
     by the ICTV. It does not compare against ALL phage genomes, just classified genomes. Having found the closet related phages 
     it runs the VIRIDIC--algorithm and parses the output to predict the taxonomy of the phage. It is only able to classify to the Genus and Species level"""
    parser = ArgumentParser(usage, description=description)
    parser.add_argument("-v", "--verbose", action="store_true", default = 0)
    parser.add_argument("-t", "--threads", dest='threads', type=str, default= "8" ,
                        help= "Maximum number of threads that will be used")
    parser.add_argument('-i', '--input', dest='in_fasta', type=str, help='Path to an input fasta file')
    parser.add_argument("-p", "--prefix",  type=str, default ="", dest='prefix',
                        help='will add the prefix to results and summary files that will store results of MASH and comparision to the VMR Data produced by'
                        'ICTV combines both sets of this data into a single csv file. '
                        'Use this flag if you want to run multiple times and keep the results files without manual renaming of files')
    parser.add_argument("-d", "--distance", type=str, default = "0.2",  dest="dist",
                        help='Will change the mash distance for the intial seraching for close relatives. We suggesting keeping at 0.2'
                        ' If this results in the phage not being classified, then increasing to 0.3 might result in an output that shows'
                        ' the phage is a new genus. We have found increasing above 0.2 does not place the query in any current genus, only'
                        ' provides the output files to demonstrate it falls outside of current genera')
    parser.add_argument("--Figures", type=str, choices=["T", "F"], default = "T",  help="Specify 'T' or 'F' to produce Figures. Using F"
                        "will speed up the time it takes to run the script - but you get no Figures. ")


    args, nargs = parser.parse_known_args()
    verbose = args.verbose

    #turn on ICECREAM reporting
    if not verbose: ic.disable()

    timer_start = time.time()

    # this is the location of where the script and the databases are (instead of current_directory which is the users current directory)
    HOME = os.path.dirname(__file__)
    VMR_path = os.path.join(HOME, 'VMR.xlsx')
    blastdb_path = os.path.join(HOME, 'Bacteriophage_genomes.fasta')
    ICTV_path = os.path.join(HOME, 'ICTV.msh')

    if os.path.exists(VMR_path):
        print_ok(f"Found {VMR_path} as expected")
    else:
        print_error(f'File {VMR_path} does not exist will try downloading now')
        print_error("Will download the current VMR now")
        url = "https://ictv.global/vmr/current"
        file_download = "VMR.xlsx"
        download_command = f"curl -o {file_download} {url}"
        try:
            subprocess.run(download_command, shell=True, check=True)
            print(f"{url} downloaded successfully!")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while downloading {url}: {e}")


    if os.path.exists(ICTV_path):
        print_ok(f" Found {ICTV_path} as expected")
    #else:
    #    print_error(f'File {ICTV_path} does not exist, was it downloaded correctly?')
    else:
        print_error(f"File {ICTV_path} does not exist will create database now  ")
        print_error("Will download the database now and create database")
        url = "https://millardlab-inphared.s3.climb.ac.uk/ICTV_2023.msh"
        file_download = "ICTV.msh"
        download_command = f"curl -o {file_download} {url}"
        try:
            subprocess.run(download_command, shell=True, check=True)
            print(f"{url} downloaded successfully!")
        except subprocess.CalledProcessError as e:
                print(f"An error occurred while downloading {url}: {e}")

    check_programs()
    check_blastDB(blastdb_path)
    
    #Defined and set some parameters
    fasta_file = args.in_fasta
    base = os.path.basename(fasta_file).removesuffix('.gz').removesuffix('.fasta').removesuffix('.fna').removesuffix('.fsa')
    cwd = os.getcwd()
    threads = args.threads
    mash_dist = args.dist
    ic(f"Number of set threads {threads}")
    #create results folder
    results_path = os.path.join(cwd, f"{base}_taxmyphage_results")
    query = os.path.join(results_path,'query.fasta')

    #path to the combined df containing mash and VMR data
    out_csv_of_taxonomy = args.prefix+"Output_of_taxonomy.csv"
    taxa_csv_output_path = os.path.join(results_path, out_csv_of_taxonomy)

    #path the final results summary file
    summary_results = args.prefix+"Summary_file.txt"
    summary_output_path = os.path.join(results_path, summary_results)

    #fasta file to store known taxa
    known_taxa_path = os.path.join(results_path, 'known_taxa.fa')
    #store files for VIRIDIC run- or equivalent
    viridic_in_path = os.path.join(results_path, 'viridic_in.fa')

    heatmap_file = os.path.join(results_path, 'heatmap')
    similarities_file = os.path.join(results_path, 'similarities.tsv')
    #Statments to output

    summary_statement1 ="""
    \n The data from the initial mash searching is below as tsv format \n
    Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing you determine if it represents a new 
    species or genus. It does not tell you if it is similar to other phages that have yet to be classified 
    You can do this by comparison with INPHARED database if you wish https://github.com/RyanCook94/inphared or BLAST etc \n\n
    """

    statement_current_genus_new_sp="""
    Query sequence can be classified within a current genus and represents a new species, it is in:\n
    """
    statement_current_genus_sp="""
    \nQuery sequence can be classified within a current genus and species, it is in:\n
    """
    summary_statement_inconsitent="""
    The number of expected genera based on current ICTV classification is less than the predicted 
    number of genus clusters as predicted by VIRIDIC-algorithm. This does not mean the current ICTV classification
    is wrong (it might be)or that VIRIDIC-algorithm is wrong. It could be an edge case that automated process cannot
    distinguish. It will require more manual curation to look at the output files
    \n 
    """


    if not os.path.exists(results_path):
        os.makedirs(results_path)
        print_ok(f"Directory '{results_path}' created to store results")
    else:
        print_warn(f"\n\nWarning: Directory '{results_path}'already exists. All results will be overwritten.")


    #write the input file to a new file with a header called "taxmyPhage" which makes it easier to find in data later on

    handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file)
    entries = list(SimpleFastaParser(handle))

    num_sequences = len(entries)

    if num_sequences == 0:
        print_error("Error: The FASTA file is empty.")
        exit()
    elif num_sequences == 1:
        # Open output FASTA file
        with open(query, "w") as output_fid:
            name, seq = entries[0]
            with open(summary_output_path, 'a') as fid:
                print (f"Query sequence header was:{name}",  file=fid)
                print (f">taxmyPhage\n{seq}", file=output_fid)
    else:
        print_error(f"\nError: The {fasta_file} FASTA file contains {num_sequences} sequences."\
              f" Only one sequence can be classified at a time ")


    #Read the viral master species record into a DataFrame
    taxa_df = pd.read_excel(VMR_path,sheet_name=0)

    #Print the DataFrame and rename a column
    ic(f"taxa_df")

    taxa_df = taxa_df.rename(columns={'Virus GENBANK accession': 'Genbank'})
    taxa_df['Genbank'].fillna('', inplace=True)
    # Get the column headings as a list

    #headings = list(taxa_df.columns)
    #create a dictionary of Accessions linking to Genus
    accession_genus_dict = taxa_df.set_index('Genbank')['Genus'].to_dict()

    #run mash to get top hit and read into a pandas dataframe
    mash_output = subprocess.check_output(['mash', 'dist', '-d' , mash_dist,  '-p', threads, ICTV_path, query])

    # list of names for the headers
    mash_df  = pd.read_csv(io.StringIO(mash_output.decode('utf-8')), sep='\t', header=None, names=['Reference', 'Query', 'distance', 'p-value', 'shared-hashes', 'ANI'])
    number_hits = mash_df.shape[0]

    #get the number of genomes wih mash distance < 0.2

    if number_hits < 1:
        print_error("""
    Error: No hits were found with the default settings
    The phage likely represents a new species and genus 
    However tax_my_phage is unable to classify it at this point as it can only classify at the Genus/Species level
              """)
        exit ()
    else:
        print_res(f"""
        Number of phage genomes detected with mash distance of < 0.2 is:{number_hits}""")

    #sort dataframe by distance so they are at the top
    mash_df = mash_df.sort_values(by='distance', ascending=True)
    mash_df.to_csv(os.path.join(results_path, 'mash.txt'), index=False)
    minimum_value = mash_df['distance'].min()
    maximum_value = mash_df.head(10)['distance'].max()

    print_ok(f"""The mash distances obtained for this query phage
    is a minimum value of {minimum_value} and maximum value of {minimum_value} """)



    #set the maximum number of hits to take forward. Max is 50 or the max number in the table if <50
    filter_hits =""
    if number_hits < 10:
        filter_hits = number_hits
    else:
        filter_hits = 10

    #copy top 50 hits to a new dataframe
    top_50 = mash_df.iloc[:filter_hits].copy()

    ic (f"{mash_df.head(10)}")
    ic (f"{top_50}")
    #reindex 
    top_50.reset_index(drop=True, inplace=True)

    value_at_50th_position = top_50['distance'].iloc[filter_hits-1]
    ic(f"{value_at_50th_position}")

    
    top_50['genus'] = top_50['Reference'].str.split('/').str[1]
    top_50['acc'] = top_50['Reference'].str.split('/').str[-1].str.split('.fna|.fsa').str[0]
    top_50 = top_50.merge(taxa_df, left_on = 'acc', right_on = 'Genbank')
    top_50['ANI'] = (1 - top_50.distance)*100

    #returns the unique genera names found in the mash hits - top_50 is not the best name!

    unique_genera_counts = top_50.genus.value_counts()
    print_ok(unique_genera_counts.to_dict())
    unique_genera = unique_genera_counts.index.tolist()

    #unique_genera top_50.genus.value_counts().to_dict()
    #print for error checking
    ic(unique_genera)

    #number of genera
    number_of_genera = len(unique_genera)

    print(f"Found {number_of_genera} genera associated with this query genome")

    #get all the keys for from a dictionary of accessions and genus names
    keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]

    #print the keys
    ic(keys)


    #Do different things depending how many unique genera were found

    if len(unique_genera) == 1:
        print_ok("Only found 1 genus so will proceed with getting all genomes associated with that genus")
        keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]
        number_ok_keys =len(keys)
        print_ok(f"Number of known species in the genus is {number_ok_keys} \n ")
        # create a command string for blastdbcmd
        get_genomes_cmd = f"blastdbcmd -db {HOME}/Bacteriophage_genomes.fasta  -entry {','.join(keys)} -out {known_taxa_path} "
        #subprocess.run(get_genomes_cmd, shell=True, check=True)
        res = getoutput(get_genomes_cmd)

    elif len(unique_genera) >1:
        print_ok("Found multiple genera that this query phage might be similar to so will proceed with processing them all")
        list_of_genus_accessions =[]
        for i in unique_genera:
            keys = [k for k, v in accession_genus_dict.items() if v == i]
            number_of_keys = len(keys)
            #ic(keys)
            list_of_genus_accessions.extend(keys)
            print_ok(f"Number of known species in the genus {i} is {number_of_keys}")
        ic(list_of_genus_accessions)
        ic(len(list_of_genus_accessions))
        get_genomes_cmd = f"blastdbcmd -db {HOME}/Bacteriophage_genomes.fasta  -entry {','.join(list_of_genus_accessions)} -out {known_taxa_path}"
        res = getoutput(get_genomes_cmd)

    #get smallest mash distance

    min_dist = top_50['distance'].min()

    if  min_dist < 0.04: 
        print_ok ("Phage is likely NOT a new species, will run further analysis now to to confirm this \n ")
        top_df = top_50[top_50['distance'] == min_dist]
        ic(top_df)


    elif min_dist > 0.04 and min_dist  < 0.1: 
        print_ok ("It is not clear if the phage is a new species or not. Will run further analysis now to confirm this...\n")
        top_df = top_50[top_50['distance'] < 0.1 ]
        ic(top_df)
        print(top_50.genus.value_counts())

    elif min_dist > 0.1  and min_dist  < 0.2:
        print_ok ("Phage is a new species. Will run further analysis now ....\n")
        top_df = top_50[top_50['distance'] < 0.1 ]
        ic(top_df)


    #######run poor mans viridic
    with open(known_taxa_path, 'r') as file1:
        with open(query, 'r') as file2:
            with open(viridic_in_path, 'w') as merged_file:
                merged_file.write(file1.read())
                merged_file.write(file2.read())

    PMV = PoorMansViridic(viridic_in_path, nthreads=threads, verbose=verbose)
    df1, pmv_outfile = PMV.run()

    # heatmap and distances
    if args.Figures != "F":
        print_ok("Will calcualte and save heatmaps now")
        heatmap(PMV.dfM, heatmap_file)
    else:
        print_error("\n Skipping calculating heatmaps and saving them \n ")

    PMV.save_similarities(similarities_file)
    
    
    taxa_df = pd.read_excel(VMR_path,sheet_name=0)

    #Print the DataFrame
    ic(taxa_df)
    #rename the columns again
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


    print(f"""\n\nTotal number of VIRIDIC-algorithm genus clusters in the input including QUERY sequence was:{total_num_viridic_genus_clusters}
    Total number of VIRIDIC-algorithm species clusters including QUERY sequence was {total_num_viridic_species_clusters} """)

    print(f"""\n\nNumber of current ICTV defined genera was :{num_unique_ICTV_genera}
    Number of VIRIDIC-algorithm predicted genera (excluding query)was :{num_unique_viridic_genus_clusters} """)


    if num_unique_ICTV_genera == num_unique_viridic_genus_clusters:
        print (f"""\n\nCurrent ICTV and VIRIDIC-algorithm predictions are consistent for the data that was used to compare against""")




    print_ok(f"Number of unique VIRIDIC-algorithm clusters at default cutoff of 70% is:{num_unique_viridic_genus_clusters}")
    print_ok(f"""Number of current ICTV genera associated with the reference genomes
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
        print ("""Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level""")

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
            print_res(f"""Query sequence is: 
                    Class:{list_of_S_data[14]}
                    Family: {list_of_S_data[15]}
                    Subfamily:{list_of_S_data[16]}
                    Genus:{list_of_S_data[17]}
                    Species:{list_of_S_data[19]}
                     """)
            with open(summary_output_path,'a') as file:
                file.write(f"""statement_current_genus_sp 
                           Class:{list_of_S_data[14]}\tFamily: {list_of_S_data[15]}\tSubfamily:{list_of_S_data[16]}\tGenus:{list_of_S_data[17]}Species:{list_of_S_data[19]}
                \n{summary_statement1 }""")

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

            print_res(f"""Query sequence is: 
            Class:{dict_exemplar_genus['Class']}
            Family: {dict_exemplar_genus['Family']} 
            Subfamily:{dict_exemplar_genus['Subfamily']}
            Genus:{dict_exemplar_genus['Genus']}
            Species:{dict_exemplar_genus['Genus']} new_name
             """)

            with open(summary_output_path, 'a') as file:
                file.write(f""" {statement_current_genus_new_sp}
    Class:{dict_exemplar_genus['Class']}\tFamily: {dict_exemplar_genus['Family']}\tSubfamily:{dict_exemplar_genus['Subfamily']}\tGenus:{dict_exemplar_genus['Genus']}\tnew_specices_name
    {summary_statement1}""")
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
        print_error(f"""{summary_statement_inconsitent}\n""")

        if query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number in list_ICTV_species_clusters:
            print_ok("""Phage is within a current genus and same as a current species 
             ....working out which one now .....""")
            if query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number in list_ICTV_species_clusters:
                print(f"""Phage is within a current genus and same as a current species 
                 ....working out which one now .....""")
                predicted_genus = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
                predicted_species = dict_species_cluster_2_species_name[query_species_cluster_number]
                print(f"""QUERY is in the genus:{predicted_genus} and is species: {predicted_species}""")
                # identify the row in the pandas data frame that is the same species
                matching_species_row = merged_df[merged_df['Species'] == predicted_species]
                ic(f"{matching_species_row}")
                list_of_S_data = matching_species_row[0:].values.flatten().tolist()
                ic(f"{list_of_S_data[14:20]}")
                ic(f"{list_of_S_data}")
                print_res(f"""Query sequence is: 
                        Class:{list_of_S_data[14]}
                        Family: {list_of_S_data[15]}
                        Subfamily:{list_of_S_data[16]}
                        Genus:{list_of_S_data[17]}
                        Species:{list_of_S_data[19]}
                         """)
                with open(summary_output_path, 'a') as file:
                    file.write(f"""statement_current_genus_sp 
                               Class:{list_of_S_data[14]}\tFamily: {list_of_S_data[15]}\tSubfamily:{list_of_S_data[16]}\tGenus:{list_of_S_data[17]}Species:{list_of_S_data[19]}
                    \n{summary_statement1}""")

                mash_df.to_csv(summary_output_path, mode='a', header=True, index=False, sep='\t')




        elif query_genus_cluster_number in list_ICTV_genus_clusters and query_species_cluster_number not in list_ICTV_species_clusters:
            print_ok("""Phage is within a current genus, BUT is representative of a new species 
                     ....working out which one now .....""")

            matching_genus_rows = merged_df[merged_df['genus_cluster'] == query_genus_cluster_number]
            dict_exemplar_genus = matching_genus_rows.iloc[0].to_dict()
            genus_value =dict_exemplar_genus['Genus']
            ic(f"f{matching_genus_rows}")
            ic(f"{genus_value}")

            print_res(f"""Query sequence is in the;
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
    {summary_statement1}
    """)
            mash_df.to_csv(summary_output_path, mode='a', header=True, index=False,sep='\t')


    run_time = str(timedelta(seconds = time.time() - timer_start))
    print(f"Run time for {fasta_file}: {run_time}", file=sys.stderr)


