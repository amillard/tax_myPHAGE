"""
This module provides functionalities for the clustering of genomes based on similarity (PMV) method in 
the taxmyphage package.
It includes the clusteringSimilarityGenome class which clusters genomes based on similarity.
"""

####################################################################################################
# Imports
####################################################################################################

import os
import glob
import gzip
import subprocess
from typing import Tuple
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from icecream import ic
# precomputed code addition
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pickle


####################################################################################################
# Classes
####################################################################################################

class ClusteringOnGenomicSimilarity:
    """
    ClusteringOnGenomicSimilarity class for clustering genomes based on similarity
    """

    def __init__(
        self,
        file: str,
        reference: str,
        # precomputed code addition
        db_dir: str,
        dont_use_precomputed: bool = False,
        genus_threshold: float = 70,
        species_threshold: float = 95,
        nthreads: int = 1,
        verbose: bool = True,
        blastn_exe: str = "blastn",
        makeblastdb_exe: str = "makeblastdb",
        output_dir: str = "",
        similarity_module: bool = False,
    ):
        """
        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class
            file (str): Path to the input fasta file
            reference (str, optional): Path to the reference fasta file. Defaults to "".
            genus_threshold (float, optional): Genus threshold. Defaults to 70.
            species_threshold (float, optional): Species threshold. Defaults to 95.
            nthreads (int, optional): Number of threads to use. Defaults to 1.
            verbose (bool, optional): Whether to print messages. Defaults to True.
            blastn_exe (str, optional): Path to the blastn executable. Defaults to "blastn".
            makeblastdb_exe (str, optional): Path to the makeblastdb executable. Defaults to "makeblastdb".
            output_dir (str, optional): Path to the output directory. Defaults to "".
            similarity_module (bool, optional): Whether to use the similarity module. Defaults to False.
            
        Returns:
            None
        """

        self.verbose = verbose
        self.file = file
        # precomputed code addition
        self.db_dir = db_dir
        self.dont_use_precomputed = dont_use_precomputed
        self.reference = reference
        self.result_dir = os.path.dirname(self.file)
        self.nthreads = nthreads
        self.genus_threshold = genus_threshold
        self.species_threshold = species_threshold
        self.blastn_exe = blastn_exe
        self.makeblastdb_exe = makeblastdb_exe
        self.pmv_outdir = output_dir
        self.pmv_outfile = ''
        self.dfT = pd.DataFrame()
        self.db_blast = reference
        self.blastn_result_file = ''
        self.similarity_module = similarity_module
        # precomputed code addition        
        self.dfM = pd.DataFrame()
        self.existing_files()

    def existing_files(self):
        """
        Checks if the input files exist

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            None
        """

        if not os.path.exists(self.file):
            raise FileNotFoundError(f"File {self.file} does not exist")

        if not os.path.exists(self.reference):
            raise FileNotFoundError(f"File {self.reference} does not exist")

    def check_if_number_of_seqs_is_more_than(self, max_n = 2):
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        handle = gzip.open(self.file, 'rt') if self.file.endswith('.gz') else open(self.file)
        for n, (name, seq) in enumerate(SimpleFastaParser(handle)):
            if n >= max_n: return True
        return False
 
    def run(self) -> Tuple[pd.DataFrame, str]:
        """
        Runs the clusteringSimilarityGenome pipeline

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            dfT (pd.DataFrame): DataFrame containing the clusters
            pmv_outfile (str): Path to the output file
        """
        ic(f"Running clusteringSimilarityGenome on {self.file}\n")
        self.makeblastdb()

        # precomputed code addition
        M_file = os.path.join(self.db_dir, 'M.pa')
        if not os.path.exists(M_file): self.dont_use_precomputed = True
        if self.similarity_module: self.dont_use_precomputed = True
        
        # if self.use_precomputed:
        #     self.use_precomputed = self.check_if_number_of_seqs_is_more_than(2)
            
        if self.dont_use_precomputed:
            print("Running full blast")
            self.blastn()
            self.parse_blastn_file()
        else:
            print("Running on pre computed")
            self.blastn_cached()
            self.parse_blastn_file()
            self.merge_cached_results()

        self.calculate_distances()
        self.cluster_all()
        return self.dfT, self.pmv_outfile

    def cluster_all(self):
        """
        Clusters all the genomes

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            None
        """

        dfTg = self.sim2cluster(self.genus_threshold, "genus")
        dfTs = self.sim2cluster(self.species_threshold, "species")
        dfT = pd.merge(dfTg, dfTs, on="genome").sort_values(
            "species_cluster genus_cluster".split()
        )
        dfT.reset_index(drop=True, inplace=True)
        self.pmv_outfile = os.path.join(
            self.result_dir, os.path.basename(self.file) + ".genus_species_clusters.tsv"
        )
        dfT.to_csv(self.pmv_outfile, index=False, sep="\t")
        self.dfT = dfT

    def sim2cluster(self, th: float, tax_level: str) -> pd.DataFrame:
        """
        Clusters the genomes based on the similarity threshold

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class
            th (float): Similarity threshold
            tax_level (str): Taxonomic level

        Returns:
            pd.DataFrame: DataFrame containing the clusters
        """

        ic("Generating graph for finding", tax_level, "clusters")
        M = self.dfM
        G = nx.from_pandas_edgelist(
            M[(M.sim >= th) & (M.A != M.B)], source="A", target="B"
        )
        singletons = list(set(M.A.unique().tolist()).difference(G.nodes()))
        G.add_nodes_from(singletons)

        graphs = [G.subgraph(x) for x in nx.connected_components(G)]
        L = []
        for n, g in enumerate(graphs):
            L.extend([(node, n + 1) for node in g.nodes()])

        return pd.DataFrame(L, columns=f"genome {tax_level}_cluster".split())

    def makeblastdb(self):
        """
        Creates the blastDB

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            None
        """

        # Find all the files created by makeblastdb and remove them
        for filename in glob.glob(f"{self.reference}*.n*"):
            os.remove(filename)

        # Get the temporary directory for tests
        if self.pmv_outdir:
            self.db_blast = os.path.join(self.pmv_outdir, "tmp_blastdb")

        cmd = f"{self.makeblastdb_exe} -in {self.reference} -out {self.db_blast} -dbtype nucl"
        ic("Creating blastn database:", cmd)
        res = subprocess.getoutput(cmd)
        ic(res)

    # precomputed code addition
    # unnecessary code duplication, could be done smarter    
    def blastn_cached(self):
        query = os.path.join(os.path.dirname(os.path.dirname(self.file)), 'query.fasta')
        outfile = os.path.join(
            self.result_dir, os.path.basename(query) + ".blastn_vs2_self.tab.gz"
        )

        if not os.path.exists(outfile):
            if self.similarity_module:
                print("ARGGGGGG, we shouldn't get here")
                cmd = (
                    f'{self.blastn_exe} -evalue 1 -max_target_seqs 10000 -num_threads {self.nthreads} '
                    f'-word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {query} '
                    f'-db {self.db_blast} '
                    '-outfmt "6 qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore" '
                    f'| sort | uniq | gzip -c > {outfile}'
                )                    
            else:
                cmd = (
                    f'{self.blastn_exe} -evalue 1 -max_target_seqs 10000 -num_threads {self.nthreads} '
                    f'-word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {query} '
                    f'-db {self.db_blast} '
                    '-outfmt "6 qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore" '
                    f'| gzip -c > {outfile}'
                )

            ic("Blasting against itself:", cmd)
            ic(cmd)
            res = subprocess.getoutput(cmd)

        self.blastn_result_file = outfile
        
    def blastn(self):
        """
        Runs blastn against itself

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            None
        """

        outfile = os.path.join(
            self.result_dir, os.path.basename(self.file) + ".blastn_vs2_self.tab.gz"
        )
        if not os.path.exists(outfile):
            if self.similarity_module:
                cmd = (
                    f'{self.blastn_exe} -evalue 1 -max_target_seqs 10000 -num_threads {self.nthreads} '
                    f'-word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {self.file} '
                    f'-db {self.db_blast} '
                    '-outfmt "6 qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore" '
                    f'| sort | uniq | gzip -c > {outfile}'
                )                    
            else:
                cmd = (
                    f'{self.blastn_exe} -evalue 1 -max_target_seqs 10000 -num_threads {self.nthreads} '
                    f'-word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {self.file} '
                    f'-db {self.db_blast} '
                    '-outfmt "6 qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore" '
                    f'| gzip -c > {outfile}'
                )

            ic("Blasting against itself:", cmd)
            ic(cmd)
            subprocess.getoutput(cmd)

        self.blastn_result_file = outfile

    def parse_blastn_file(self):
        """
        Parses the blastn file

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            None
        """

        ic("Reading", self.blastn_result_file)

        num_lines = rawgencount(self.blastn_result_file)

        self.size_dict = {}
        M = {}

        previous_pair = ""

        with gzip.open(self.blastn_result_file, "rt") as df:
            genome_name = os.path.dirname(self.file).split("/")[-1]

            for line in tqdm(
                df, desc=f"{genome_name}: Blast reading:", total=num_lines, leave=False
            ):
                # do something with the line
                (
                    qseqid,
                    sseqid,
                    pident,
                    length,
                    qlen,
                    slen,
                    mismatch,
                    nident,
                    gapopen,
                    qstart,
                    qend,
                    sstart,
                    send,
                    qseq,
                    sseq,
                    evalue,
                    bitscore,
                ) = line.rstrip().split()
                key = (qseqid, sseqid)

                # if the key is different from the previous one, convert identity vector to identity values
                if key != previous_pair:
                    # only do it if the previous key is not empty (first iteration)
                    if previous_pair:
                        M[previous_pair] = np.where(M[previous_pair] != 0, 1, 0)
                        M[previous_pair] = np.sum(M[previous_pair])

                    previous_pair = key

                M.setdefault(key, np.zeros(int(qlen)))

                if qseqid not in self.size_dict:
                    self.size_dict[qseqid] = int(qlen)
                if sseqid not in self.size_dict:
                    self.size_dict[sseqid] = int(slen)

                # convert the strings to numpy arrays
                qseq = np.frombuffer(qseq.encode("utf-8"), dtype="S1")
                sseq = np.frombuffer(sseq.encode("utf-8"), dtype="S1")

                v = np.where(qseq == sseq, 1, 0)

                # find the indices of elements that are not equal to '-'. Here it is b'-' because the array is of type bytes.
                idx = qseq != b"-"

                # add the values to the matrix
                M[key][int(qstart) - 1 : int(qend)] += v[idx]

        # precomputed code addition
        # I think this is needed also in the non-precomputed
        if previous_pair:
            # Convert the last pair of the matrix to identity values
            M[previous_pair] = np.where(M[previous_pair] != 0, 1, 0)
            M[previous_pair] = np.sum(M[previous_pair])

        self.M = M

    # precomputed code addition
    def merge_cached_results(self):
        M_file = os.path.join(self.db_dir, 'M.pa')
        
        file = self.reference
        handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
        entries = list(SimpleFastaParser(handle))

        self.size_dict = {name:len(seq) for name, seq in entries}
        names = list(self.size_dict.keys())
        dfM = pd.read_parquet(M_file)
        
            
        genome_arr = np.array(list(self.M.keys()))
        dfMs = pd.DataFrame(genome_arr, columns=["A", "B"])
        dfMs["idAB"] = self.M.values()
        dfMs = dfMs['A B idAB'.split()]
        dfM = pd.concat([dfM['A B idAB'.split()], dfMs['A B idAB'.split()]])
        dfM = dfM[dfM.A.isin(names) & dfM.B.isin(names)]

        self.dfM = dfM

    def calculate_distances(self):
        """
        Calculates the distances between genomes

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            None
        """
        # precomputed code addition
        if self.dfM.empty:
            M = self.M
            genome_arr = np.array(list(M.keys()))

            dfM = pd.DataFrame(genome_arr, columns=["A", "B"])

            dfM["idAB"] = M.values()
        else:
            dfM = self.dfM
        
        M = self.M
        size_dict = self.size_dict


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
        dfM["lA"] = dfM["A"].map(size_dict)
        dfM["lB"] = dfM["B"].map(size_dict)

        # Calculate the similarity
        dfM["simAB"] = ((dfM.idAB + dfM.idBA) * 100) / (dfM.lA + dfM.lB)

        # Calculate the distance
        dfM["distAB"] = 100 - dfM.simAB

        # Calculate the aligned fraction of the genome
        dfM["afg1"] = dfM.idAB / dfM.lA
        dfM["afg2"] = dfM.idBA / dfM.lB
        dfM["glr"] = dfM[["lA", "lB"]].min(axis=1) / dfM[["lA", "lB"]].max(axis=1)

        # Calculate the similarity
        dfM["sim"] = 100 - dfM.distAB

        # Remove the duplicate pairs
        dfM["ordered_pair"] = dfM.apply(lambda x: str(sorted(x.pair_BA)), axis=1)
        dfM = dfM.drop_duplicates("ordered_pair").reset_index(drop=True)

        # Remove the columns that are not needed
        dfM = dfM.drop(
            columns=[
                "pair_BA",
                "idAB",
                "idBA",
                "lA",
                "lB",
                "simAB",
                "ordered_pair",
            ]
        )

        self.dfM = dfM

    def save_similarities(self, outfile: str = "similarities.tsv"):
        """
        Saves the similarities

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class
            outfile (str, optional): Path to the output file. Defaults to "similarities.tsv".

        Returns:
            None
        """

        df = self.dfM[["A", "B", "sim"]]
        df = df[df.A != df.B]
        df.sort_values("sim", ascending=False, inplace=True)

        # Save the dataframe
        df.to_csv(outfile, index=False, sep="\t")

    def get_query_closest(self):
        """
        Gets the query closest genomes

        Args:
            self (clusteringSimilarityGenome): clusteringSimilarityGenome class

        Returns:
            pd.DataFrame: DataFrame containing the query similarities
        """

        df = self.dfM[["A", "B", "sim"]]
        df = df[df.A != df.B]
        df.sort_values("sim", ascending=False, inplace=True)

        df_query = df[(df.A.str.contains("query_")) | (df.B.str.contains("query_"))].reset_index(drop=True)

        closest_genomes = df_query.A.values[0] if not df_query.A.values[0].startswith("query_") else df_query.B.values[0]

        closest_genomes = closest_genomes.split(":")[0]

        return closest_genomes

####################################################################################################
# Functions
####################################################################################################


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
        b = reader(1024 * 1024)


####################################################################################################


def rawgencount(filename: str) -> int:
    """Count the number of lines in a file.
    Args:
        filename (str): The name of the file to count.
    Returns:
        int: The number of lines in the file.
    """
    f = gzip.open(filename, "rb")
    f_gen = _make_gen(f.read)
    return sum(buf.count(b"\n") for buf in f_gen)


####################################################################################################
