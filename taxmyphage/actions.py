"""
This module provides functionalities for classifying bacteriophages based on their genomes.
It includes functions for running all classification methods, generating heatmaps, and handling 
files.
"""

####################################################################################################
# Imports
####################################################################################################


import os
import sys
import time
import io
from datetime import timedelta
from argparse import Namespace
from icecream import ic
from tqdm import tqdm
from Bio import SeqIO
import pandas as pd

from taxmyphage.utils import print_error, print_ok, create_folder, print_warn
from taxmyphage.pmv import ClusteringOnGenomicSimilarity
from taxmyphage.plot import heatmap
from taxmyphage.handle_files import create_files_and_result_paths, read_write_fasta
from taxmyphage.classify import (
    classification_mash,
    classification_similarity,
    classification,
)


####################################################################################################
# Functions
####################################################################################################


def all_classification(
    args: Namespace,
    taxa_df: pd.DataFrame,
    mash_dist: int,
    mash_index_path: str,
    blastdb_path: str,
    threads: int,
    verbose: bool,
    max_entries: int=0,
    force_overwrite: bool=0,
        ) -> pd.DataFrame:
    """
    This function will run all the classification methods and return the result

    Args:
        args (Namespace): arguments from the command line
        taxa_df (DataFrame): dataframe containing the taxonomy of the VMR
        mash_dist (int): distance to use for the mash analysis
        mash_index_path (str): path to the mash index
        blastdb_path (str): path to the blastdb
        threads (int): number of threads to use
        verbose (bool): if True, print more information
        max_entries (int): max number of sequences to analyse (default is no limit)
    Returns:
        dataframe with taxonomy results
    """

    tmp_fasta = os.path.join(args.output, "tmp.fasta")
    # Create a multifasta file to parse line by line

    # create a dictionary to store the taxonomy of each genome
    dict_taxonomy = {}

    num_genomes = create_files_and_result_paths(args.in_fasta, tmp_fasta)
    if max_entries and num_genomes > max_entries:
        print_error(f"{num_genomes} entries in fasta file larger than allowed number ({max_entries}. Skipping {num_genomes-max_entries}")
    else:
        print("\n------------------------------")
        print("Starting tax_my_phage analysis")
        print("------------------------------\n")

    parser = SeqIO.parse(tmp_fasta, "fasta")

    genome_counter = 0
    for genome in tqdm(parser, desc="Classifying", total=num_genomes):
        genome_counter +=1
        if max_entries and genome_counter > max_entries: break
        genome_id = genome.id

        # replace any characters that are not allowed in a file name
        for char in [
            " ",
            "/",
            "|",
            ":",
            "(",
            ")",
            "[",
            "]",
            "{",
            "}",
            "<",
            ">",
            "#",
            "%",
            "&",
            "+",
            "$",
            "=",
        ]:
            if char in genome_id:
                genome_id = genome_id.replace(char, "_")

        results_path = os.path.join(args.output, "Results_per_genome", genome_id)

        if os.path.exists(results_path):
            print(f"\n\nResult folder {results_path} exists.\n")
            if force_overwrite:
                print(f"\n\nForcefully overwriting {results_path}.\n")
            else:
                print(f"\n\nTo overwrite, add the -f flag. Ignoring {genome_id}.\n")
                continue
            
        print(f"\n\nClassifying {genome.id} in result folder {results_path}\n")

        timer_start = time.time()


        # create the results folder
        create_folder(results_path)

        # create the path to the query fasta file
        query = os.path.join(results_path, "query.fasta")

        # create a fasta file with just the query genome and add query_ to the id
        with open(query, "w", encoding="utf-8") as output_fid:
            genome.name = genome.description = ""
            genome.id = f"query_{genome_id}"
            SeqIO.write(genome, output_fid, "fasta")

        # path to the combined df containing mash and VMR data
        out_csv_of_taxonomy = args.prefix + "Output_of_taxonomy.tsv"
        taxa_csv_output_path = os.path.join(results_path, out_csv_of_taxonomy)

        # fasta file to store known taxa
        known_taxa_path = os.path.join(results_path, "known_taxa.fa")

        print("-------------")
        print("MASH analysis")
        print("-------------\n")

        mash_df, accession_genus_dict = classification_mash(
            known_taxa_path=known_taxa_path,
            results_path=results_path,
            dist=mash_dist,
            query=query,
            mash_index_path=mash_index_path,
            blastdb_path=blastdb_path,
            taxa_df=taxa_df,
            taxa_csv_output_path=taxa_csv_output_path,
            threads=threads,
            mash_exe=args.mash,
            blastdbcmd_exe=args.blastdbcmd,
        )

        if mash_df.empty:
            dict_taxonomy[genome_id] = {
                "Realm": "Unknown",
                "Kingdom": "Unknown",
                "Phylum": "Unknown",
                "Class": "Unknown",
                "Order": "Unknown",
                "Family": "Unknown",
                "Subfamily": "Unknown",
                "Genus": "New_genus",
                "Species": "New_species",
                "Message": "No hits were found with the default settings",
            }
            continue

        print("-----------------------------------------")
        print("Clustering on genomic similarity analysis")
        print("-----------------------------------------")

        merged_df, query_merged_df, closest_genome = classification_similarity(
            known_taxa_path=known_taxa_path,
            query=query,
            # precomputed code addition
            db_dir=os.path.dirname(blastdb_path),
            dont_use_precomputed=args.no_precomputed,
            taxa_df=taxa_df,
            taxa_csv_output_path=taxa_csv_output_path,
            results_path=results_path,
            threads=threads,
            accession_genus_dict=accession_genus_dict,
            Figure=args.Figure,
            verbose=verbose,
            blastn_exe=args.blastn,
            makeblastdb_exe=args.makeblastdb,
        )

        print("--------------------")
        print("Final classification")
        print("--------------------")

        genome_taxo = classification(
            merged_df=merged_df,
            query_merged_df=query_merged_df,
            results_path=results_path,
            mash_df=mash_df,
            prefix=args.prefix,
            closest_genome=closest_genome,
        )

        dict_taxonomy[genome_id] = genome_taxo

        run_time = str(timedelta(seconds=time.time() - timer_start))
        print(f"Run time for {genome_id}: {run_time}\n", file=sys.stderr)
        print("-" * 80, file=sys.stderr)

    # write the taxonomy to a csv file
    taxonomy_tsv = os.path.join(args.output, "Summary_taxonomy.tsv")

    txt = "Genome\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tSubfamily\tGenus\tSpecies\tFull_taxonomy\tMessage\n"

    for key, value in dict_taxonomy.items():
        string_taxo = ""
        full_string = ""
        for taxo in [
            "Realm",
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Subfamily",
            "Genus",
            "Species",
        ]:
            # Change the taxonomy to Not Defined Yet if it is empty or nan
            taxonomy = (
                value[taxo]
                if value[taxo] != "" and value[taxo] == value[taxo]
                else "Not Defined Yet"
            )
            string_taxo += f"{taxonomy}\t"

            # Change the taxonomy to empty string if it is empty or nan
            taxonomy_full = (
                value[taxo]
                if value[taxo] != "" and value[taxo] == value[taxo]
                else ""
            )
            prefix = "sf" if taxo == "Subfamily" else taxo[0].lower()
            full_string += f"{prefix}__{taxonomy_full};"

        # remove the last tab
        string_taxo = string_taxo.rstrip("\t")
        full_string = full_string[:-1]

        # Remove query_ from the genome id
        query = key.replace("query_", "")

        # Get the information message
        info = value["Message"]

        txt = txt + f"{query}\t{string_taxo}\t{full_string}\t{info}\n"

    dfT = pd.read_table(io.StringIO(txt))
    dfT.to_csv(taxonomy_tsv, sep='\t')

    # clean up
    os.remove(tmp_fasta)
    return dfT


####################################################################################################


def clustering_on_genomic_similarity(args: Namespace, threads: int, verbose: bool) -> None:
    """
    Run the clustering on genomic similarity algorithm on the input multi fasta file

    Args:
        args (Namespace): arguments from the command line
        threads (int): number of threads to use
        verbose (bool): if True, print more information

    Returns:
        None
    """

    print("Similarity analysis...")

    tmp_fasta = os.path.join(args.output, "pmv.fasta")

    num_genomes = create_files_and_result_paths(args.in_fasta, tmp_fasta)
    if num_genomes <= 1:
        print_error("Cannot create similarity matrix with one genome only.")
        return
    
    # Check if reference is provided
    if args.reference:
        print_ok(f"\nUsing {args.reference} as reference\n")

        reference = os.path.join(args.output, "reference_pmv.fasta")

        with open(reference, "wt", encoding="utf-8") as f:
            read_write_fasta(args.reference, f)

        # Not possible if matrix not square at the moment
        args.Figure = False
    else:
        print_warn(
            "No reference provided. Will use the first genome in the input file as reference"
        )
        parser = SeqIO.parse(tmp_fasta, "fasta")
        reference = tmp_fasta

    results_path = os.path.join(args.output)

    heatmap_file = os.path.join(results_path, "heatmap")
    top_right_matrix = os.path.join(results_path, "top_right_matrix.tsv")
    similarities_file = os.path.join(results_path, "similarities.tsv")

    print_ok(f"\nCalculating the genomic similarity in {results_path}...")

    # run ClusteringOnGenomicSimilarity
    pmv = ClusteringOnGenomicSimilarity(
        file=tmp_fasta,
        reference=reference,
        db_dir='',
        dont_use_precomputed=False,
        nthreads=threads,
        verbose=verbose,
        blastn_exe=args.blastn,
        makeblastdb_exe=args.makeblastdb,
        similarity_module=True,
    )

    dfT, pmv_outfile = pmv.run()

    # heatmap and distances
    if args.Figure:
        print_ok("\nWill calculate and save heatmaps now\n")
        accession_genus_dict = {name: "" for name in pmv.dfM.A.unique()}
        heatmap(pmv.dfM, heatmap_file, top_right_matrix, accession_genus_dict)
    else:
        print_error("\nSkipping calculating heatmaps and saving them\n ")

    pmv.save_similarities(similarities_file)

    df = pmv.dfM
    df = df[df.A != df.B].set_index(["A", "B"])

    # print(df.to_csv(sep='\t', float_format='%.4f'))
    ic(df)
