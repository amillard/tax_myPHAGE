#!/usr/bin/env python3

'''
Accessory script
Using Thomas pmv class
Allows analysis of large number of genomes and clustering into genus/species groups
'''

import subprocess
import sys
import os
import io
import gzip
import time
import re
import glob
from argparse import ArgumentParser
from itertools import combinations_with_replacement
import numpy as np
import pandas as pd
from icecream import ic
from Bio.SeqIO.FastaIO import SimpleFastaParser
import networkx as nx
from tqdm import tqdm
from datetime import timedelta


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
        dfTg = self.sim2cluster(self.genus_threshold, 'genus')
        dfTs = self.sim2cluster(self.species_threshold, 'species')
        dfT = pd.merge(dfTg, dfTs, on='genome').sort_values('species_cluster genus_cluster'.split())
        dfT.reset_index(drop=True, inplace=True)
        self.pmv_outfile = os.path.join(self.result_dir, os.path.basename(self.file) + '.pmv_genus_species_clusters.tsv')
        dfT.to_csv(self.pmv_outfile, index=False, sep='\t')
        self.dfT = dfT

    def sim2cluster(self, th, tax_level):
        ic("Generating graph for finding", tax_level, "clusters")
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
            res = subprocess.getoutput(cmd)
        
    def blastn(self):
        outfile = os.path.join(self.result_dir, os.path.basename(self.file) + '.blastn_vs2_self.tab.gz')
        if not os.path.exists(outfile):
            cmd = f'blastn -evalue 1 -max_target_seqs 10000 -num_threads {self.nthreads} -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {self.file} -db {self.file} -outfmt "6 qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore" | gzip -c > {outfile}'
            ic("Blasting against itself:", cmd)
            subprocess.getoutput(cmd)

        self.blastn_result_file = outfile

    def parse_blastn_file(self):
        # 10 times faster parsing function from Remi Denise
        ic("Reading", self.blastn_result_file)
        
        # helper functions for tqdm progress bar
        def _make_gen(reader):
            """Generator to read a file piece by piece.
            Default chunk size: 1k.
            Args:
                reader (func): Function to read a piece of the file.
            Yields:
                generator: A generator object that yields pieces of the file.
            """
            b = reader(1024 * 1024)
            while b:
                yield b
                b = reader(1024*1024)

        def rawgencount(filename):
            """Count the number of lines in a file.
            Args:
                filename (str): The name of the file to count.
            Returns:
                int: The number of lines in the file.
            """
            f = gzip.open(filename, 'rb')
            f_gen = _make_gen(f.read)
            return sum(buf.count(b'\n') for buf in f_gen)
        
        num_lines = rawgencount(self.blastn_result_file)

        self.size_dict = {}
        M = {}
        with gzip.open(self.blastn_result_file, 'rt') as df:
            for line in tqdm(df, total=num_lines):
                # do something with the line
                qseqid, sseqid, pident, length, qlen, slen, mismatch, nident, gapopen, qstart, qend, sstart, send, qseq, sseq, evalue, bitscore = line.rstrip().split()
                key = (qseqid, sseqid)
                M.setdefault(key, np.zeros(int(qlen)))

                if qseqid not in self.size_dict:
                    self.size_dict[qseqid] = int(qlen) 
                if sseqid not in self.size_dict:
                    self.size_dict[sseqid] = int(slen)

                # convert the strings to numpy arrays
                qseq = np.frombuffer(qseq.encode('utf-8'), dtype="S1")
                sseq = np.frombuffer(sseq.encode('utf-8'), dtype="S1")

                v = np.where(qseq == sseq, 1, 0)

                # find the indices of elements that are not equal to '-'. Here it is b'-' because the array is of type bytes.
                idx = qseq != b'-'

                # add the values to the matrix
                M[key][int(qstart)-1:int(qend)] += v[idx]

        # convert the matrix to a binary matrix
        M = {key:np.where(value != 0, 1, 0) for key, value in M.items()}

        self.M = M
    
    def calculate_distances(self):
        M = self.M
        size_dict = self.size_dict

        genome_arr = np.array(list(M.keys()))
        identity_arr = np.array(list(M.values()), dtype=object).reshape(-1, 1)
        dfM = pd.DataFrame(np.hstack([genome_arr, identity_arr]), columns=['A', 'B', 'identity_seq'])

        dfM["idAB"] = dfM['identity_seq'].apply(lambda x: np.sum(x))

        # creating a dictionary of genome name identity 
        # As the blast is double sided need to check the identity of both genomes by looking at the opposite pair
        dict_BA = dfM.set_index(["A", "B"]).idAB.to_dict()

        # Creating the pair of genomes in order B, A
        dfM["pair_BA"] = dfM.apply(lambda x: (x.B, x.A), axis=1)

        # Setting the identity of the pair B, A
        dfM["idBA"] = dfM.pair_BA.map(dict_BA)

        # If the identity of the pair B, A is NaN then the pair is A, B
        dfM.loc[dfM.idBA.isna(), "idBA"] = dfM.loc[dfM.idBA.isna(), "idAB"]

        # Map the size of the genome to the dataframe
        dfM["lA"] = dfM['A'].map(size_dict)
        dfM["lB"] = dfM['B'].map(size_dict)

        # Calculate the similarity
        dfM["simAB"] = ((dfM.idAB + dfM.idBA) * 100) / (dfM.lA + dfM.lB)

        # Calculate the distance
        dfM["distAB"] = 100 - dfM.simAB

        # Calculate the aligned fraction of the genome
        dfM["afg1"] = dfM.idAB / dfM.lA
        dfM["afg2"] = dfM.idBA / dfM.lB
        dfM["glr"] = dfM[["lA", "lB"]].min(axis=1) / dfM[["lA", "lB"]].max(axis=1)

        # Calculate the similarity
        dfM['sim'] = 100 - dfM.distAB

        # Remove the duplicate pairs
        dfM["ordered_pair"] = dfM.apply(lambda x: str(sorted(x.pair_BA)), axis=1)
        dfM = dfM.drop_duplicates("ordered_pair").reset_index(drop=True)

        # Remove the columns that are not needed
        dfM = dfM.drop(columns=['identity_seq', 'pair_BA', 'idAB', 'idBA', 'lA', 'lB', 'simAB', 'ordered_pair'])

        self.dfM = dfM
    
    def save_similarities(self, outfile='pmv_similarities.tsv'):
        df = self.dfM['A B sim'.split()]
        df = df[df.A != df.B]
        df.sort_values('sim', ascending=False, inplace=True)
        df.index.name = ''
        df.to_csv(outfile, index=False, sep='\t')


if __name__ == '__main__': 
usage = "%prog [options] file (or - for stdin)"
description= """Provide a multi fasta input file to run a VIRIDIC like clustering on it
Does not link to taxonomy
Provides an output file of cluster numbers for Genus & species """
parser = ArgumentParser(usage, description=description)
parser.add_argument("-t", "--threads", dest='threads', type=str, default= "8",
                    help= "Maximum number of threads that will be used")
parser.add_argument('-i', '--input', dest='in_fasta', type=str, help='Path to fasta file')

parser.add_argument('-o', "--outdir", type=str, default = None, help="Change the path to the output directory")
args = parser.parse_args()
print (f"{args.in_fasta}")

mypmv = PoorMansViridic(args.in_fasta)

mypmv.run()

mypmv.save_similarities()
