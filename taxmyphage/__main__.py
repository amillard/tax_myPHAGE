#!/usr/bin/env python3

"""
This module is used for classifying bacteriophages based on their genomes.
It provides functionalities for checking and installing necessary databases,
handling files, and performing classifications.
"""

###################################################################################################
# Imports
###################################################################################################


import os
from icecream import ic

from taxmyphage import cli
from taxmyphage.download_check import (
    check_blastDB,
    check_mash_index,
    check_blastn_dataframe,
    check_VMR,
    install_db,
)

from taxmyphage import __version__
from taxmyphage.utils import create_folder, CheckSoftware
from taxmyphage.handle_files import read_VMR
from taxmyphage.actions import all_classification, clustering_on_genomic_similarity


###################################################################################################
# Functions
###################################################################################################


def main():
    """
    This is the main function of the script
    """
    # Set up the arguments
    args, nargs = cli.cli()

    verbose = args.verbose

    print("="*25)
    print(f"TAXMYPHAGE version {__version__}")
    print("="*25 + "\n")

    # turn on ICECREAM reporting
    if not verbose:
        ic.disable()

    # this is the location of where the script and the databases are stored
    vmr_path = os.path.join(args.db_folder, "VMR.xlsx")
    blastdb_path = os.path.join(args.db_folder, "Bacteriophage_genomes.fasta")
    mash_index_path = os.path.join(args.db_folder, "ICTV_2023.msh")
    blastn_df_path = os.path.join(args.db_folder, "M.pa")
    
    # Check if the user wants to install the database
    if args.command == "install":
        CheckSoftware(args.makeblastdb)

        install_db(
            VMR_path=vmr_path,
            blastdb_path=blastdb_path,
            mash_index_path=mash_index_path,
            blastn_df_path=blastn_df_path,
            output=args.db_folder,
            makeblastdb=args.makeblastdb,
        )

    list_softwares = [
        "args.blastdbcmd",
        "args.blastn",
        "args.makeblastdb",
        "args.mash",
    ]

    list_softwares_names = [
        "blastdbcmd",
        "blastn",
        "makeblastdb",
        "mash",
    ]

    software2check = [eval(list_softwares[software]) for software in range(4) if list_softwares_names[software] in args]

    for exe in software2check:
        CheckSoftware(exe)

    # Defined and set some parameters
    threads = args.threads

    create_folder(args.output)

    if args.command == "run":
        print("--------------------------")
        print("Looking for database files")
        print("--------------------------\n")

        # Get the mash distance
        mash_dist = args.dist

        # Check if the VMR file exists
        check_VMR(vmr_path)

        # read in the VMR data
        taxa_df = read_VMR(VMR_path=vmr_path)

        # Check if the mash index exists
        check_mash_index(mash_index_path)

        # Check if the precomputed blastn dataframe exists
        check_blastn_dataframe(blastn_df_path)

        # Check if the blastDB file exists
        all_phages_name = check_blastDB(
            blastdb_path,
            output=args.output,
            makeblastdb_exe=args.makeblastdb,
        )

        # Reduce VMR to only the phages in the blastDB
        taxa_df = taxa_df[taxa_df["Genbank"].isin(all_phages_name)].reset_index(
            drop=True
        )

        # Run the classification
        all_classification(
            args=args,
            taxa_df=taxa_df,
            threads=threads,
            blastdb_path=blastdb_path,
            mash_index_path=mash_index_path,
            mash_dist=mash_dist,
            verbose=verbose,
            force_overwrite=args.force_overwrite,
        )

    elif args.command == "similarity":
        # Run Clustering on Genomic Similarity
        clustering_on_genomic_similarity(
            args=args,
            threads=threads,
            verbose=verbose,
        )


###################################################################################################

if __name__ == "__main__":
    main()

###################################################################################################
