#!/usr/bin/env python3

############################################################################################################
# Imports
############################################################################################################


import os
from icecream import ic
from tqdm import tqdm
from Bio import SeqIO
import time

# Import the PoorMansViridic class
from taxmyphage import cli
from taxmyphage.download_check import check_blastDB, check_mash_index, check_VMR
from taxmyphage.utils import create_folder, print_ok
from taxmyphage.handle_files import create_files_and_result_paths, read_VMR
from taxmyphage.classify import (
    classification_mash,
    classification_viridic,
    classification,
)

############################################################################################################

if __name__ == "__main__":
    # Set up the arguments
    args, nargs = cli.cli()

    verbose = args.verbose
    # Defined and set some parameters
    threads = args.threads

    mash_dist = args.dist

    create_folder(args.output)

    # turn on ICECREAM reporting
    if not verbose:
        ic.disable()

    # this is the location of where the script and the databases are (instead of current_directory which is the users current directory)
    VMR_path = args.VMR_file
    blastdb_path = args.ICTV_db
    mash_index_path = args.mash_index

    print("Looking for database files...\n")

    # Check if the VMR file exists
    check_VMR(VMR_path)

    # read in the VMR data
    taxa_df = read_VMR(VMR_path=VMR_path)

    # Check if the mash index exists
    check_mash_index(mash_index_path)

    # Check if the blastDB file exists
    check_blastDB(blastdb_path, output=args.output)

    tmp_fasta = os.path.join(args.output, "tmp.fasta")
    # Create a multifasta file to parse line by line
    num_genomes = create_files_and_result_paths(args.in_fasta, tmp_fasta)

    parser = SeqIO.parse(tmp_fasta, "fasta")

    for genome in tqdm(parser, desc="Classifying", total=num_genomes):
        results_path = os.path.join(args.output, genome.id)
        print_ok(f"\nClassifying {genome.id} in result folder {results_path}...")

        timer_start = time.time()

        print("\nStarting tax_my_phage analysis...\n")

        # create the results folder
        create_folder(results_path)

        # create the path to the query fasta file
        query = os.path.join(results_path, "query.fasta")

        # create a fasta file with just the query genome and add query_ to the id
        with open(query, "w") as output_fid:
            genome.name = genome.description = ""
            genome.id = f"query_{genome.id}"
            SeqIO.write(genome, output_fid, "fasta")

        # path to the combined df containing mash and VMR data
        out_csv_of_taxonomy = args.prefix + "Output_of_taxonomy.csv"
        taxa_csv_output_path = os.path.join(results_path, out_csv_of_taxonomy)

        # fasta file to store known taxa
        known_taxa_path = os.path.join(results_path, "known_taxa.fa")

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
        )

        merged_df, copy_merged_df = classification_viridic(
            known_taxa_path=known_taxa_path,
            query=query,
            taxa_df=taxa_df,
            taxa_csv_output_path=taxa_csv_output_path,
            results_path=results_path,
            threads=threads,
            accession_genus_dict=accession_genus_dict,
            Figure=args.Figure,
            verbose=verbose,
        )

        classification(
            merged_df=merged_df,
            copy_merged_df=copy_merged_df,
            results_path=results_path,
            mash_df=mash_df,
            prefix=args.prefix,
            genome=genome,
            timer_start=timer_start,
        )

    # clean up
    os.remove(tmp_fasta)
