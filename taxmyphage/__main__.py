#!/usr/bin/env python3

############################################################################################################
# Imports
############################################################################################################


import os, sys
from icecream import ic
from tqdm import tqdm
from Bio import SeqIO
import time
from datetime import timedelta

from taxmyphage import cli
from taxmyphage.download_check import (
    check_blastDB,
    check_mash_index,
    check_VMR,
    install_db,
)
from taxmyphage.utils import create_folder, print_ok, CheckSoftware
from taxmyphage.handle_files import create_files_and_result_paths, read_VMR
from taxmyphage.actions import all_classification, viridic
from taxmyphage.classify import (
    classification_mash,
    classification_viridic,
    classification,
)

############################################################################################################
# Functions
############################################################################################################


def main():
    # Set up the arguments
    args, nargs = cli.cli()

    verbose = args.verbose

    # turn on ICECREAM reporting
    if not verbose:
        ic.disable()

    # this is the location of where the script and the databases are stored
    VMR_path = os.path.join(args.db_folder, "VMR.xlsx")
    blastdb_path = os.path.join(args.db_folder, "Bacteriophage_genomes.fasta")
    mash_index_path = os.path.join(args.db_folder, "ICTV_2023.msh")

    # Check if the user wants to install the database
    if args.command == "install":
        CheckSoftware(args.makeblastdb)

        install_db(
            VMR_path=VMR_path,
            blastdb_path=blastdb_path,
            mash_index_path=mash_index_path,
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
        print("Looking for database files...\n")

        # Get the mash distance
        mash_dist = args.dist

        # Check if the VMR file exists
        check_VMR(VMR_path)

        # read in the VMR data
        taxa_df = read_VMR(VMR_path=VMR_path)

        # Check if the mash index exists
        check_mash_index(mash_index_path)

        # Check if the blastDB file exists
        check_blastDB(
            blastdb_path,
            output=args.output,
            makeblastdb_exe=args.makeblastdb,
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
        )

    elif args.command == "viridic":
        # Run Poor Man VIRIDIC
        viridic(
            args=args,
            threads=threads,
            verbose=verbose,
        )


############################################################################################################

if __name__ == "__main__":
    main()

############################################################################################################
