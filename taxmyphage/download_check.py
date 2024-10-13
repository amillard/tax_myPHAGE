"""
This module provides functionalities for checking and downloading necessary databases for the 
taxmyphage package.
It includes functions to check if the necessary databases exist and to download them if they do not.
"""

####################################################################################################
# Imports
####################################################################################################

import os
import sys
import subprocess
import gzip
import shutil
from typing import List
import wget
from Bio import SeqIO

from taxmyphage.utils import print_error, print_ok, create_folder

####################################################################################################
# Functions
####################################################################################################


def download(url: str, output: str) -> None:
    """
    Downloads a file from a URL to a specified output path

    Args:
        url (str): URL to download from
        output (str): Path to the output file

    Returns:
        None
    """

    try:
        wget.download(url, output)
        print(f"\n{url} downloaded successfully!")
    except Exception as e:
        print(f"An error occurred while downloading {url}: {e}")
        sys.exit()


####################################################################################################


def makeblastdb(blastdb_path: str, makeblastdb_exe: str) -> None:
    """
    Creates the blastDB

    Args:
        blastdb_path (str): Path to the blastDB
        makeblastdb_exe (str): Path to the makeblastdb executable

    Returns:
        None
    """

    makeblastdb_command = (
        f"{makeblastdb_exe} -in {blastdb_path} -parse_seqids -dbtype nucl"
    )
    try:
        subprocess.run(makeblastdb_command, shell=True, check=True)
        print("makeblastdb command executed successfully!\n")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing makeblastdb: {e}")
        sys.exit()


####################################################################################################


def unzip_file(file_path: str, output_path: str) -> None:
    """
    Unzips a file

    Args:
        file_path (str): Path to the file to unzip
        output_path (str): Path to the output file

    Returns:
        None
    """

    try:
        with gzip.open(f"{file_path}", "rb") as f_in:
            with open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        print("File unzipped successfully!")
    except Exception as e:
        print(f"An error occurred while unzipping the file: {e}")
        sys.exit()


####################################################################################################


def check_blastDB(
    blastdb_path: str, output: str, makeblastdb_exe: str, install: bool = False
) -> List[str]:
    """
    Checks if the blastDB is present and if not downloads it and creates the database

    Args:
        blastdb_path (str): Path to the blastDB
        output (str): Path to the output directory
        makeblastdb_exe (str): Path to the makeblastdb executable
        install (bool, optional): Whether to install the blastDB. Defaults to False.

    Returns:
        List[str]: List of all phage names in the blastDB
    """

    # Variable to check if the uncompressed blastDB needs to be removed to save space
    to_remove = False

    # Check if the blastDB has been created and if not create it
    if os.path.exists(blastdb_path + ".nhr"):
        print_ok(f"Found {blastdb_path}.nhr as expected")
    else:
        # check if blastDB is present
        if os.path.exists(blastdb_path):
            print_ok(f"Found {blastdb_path} as expected")

            # Unzip the file if it is zipped
            if blastdb_path.endswith(".gz"):
                # Get the basename of the blastdb_path without the .gz extension
                blastdb_path_no_gz = os.path.basename(blastdb_path[:-3])
                blastdb_path_no_gz = os.path.join(output, blastdb_path_no_gz)

                unzip_file(blastdb_path, blastdb_path_no_gz)
                blastdb_path = blastdb_path_no_gz
                to_remove = True
        else:
            if install:
                print_error(
                    f"File {blastdb_path} does not exist will create database now  "
                )
                print_error("Will download the database now and create database")

                url = "https://millardlab-taxmyphage.s3.climb.ac.uk/Bacteriophage_genomes_MSL39v1.fasta.gz"

                create_folder(os.path.dirname(blastdb_path))

                # Download the file from the URL to the output directory
                download(url, f"{blastdb_path}.gz")

                # Gunzip the file
                unzip_file(f"{blastdb_path}.gz", blastdb_path)

                to_remove = True
            else:
                print_error(
                    f"File {blastdb_path} does not exist. Please download the database and create the database."
                    "Or use the --install flag to download and create the database."
                )
                sys.exit()

        makeblastdb(blastdb_path, makeblastdb_exe)

        if to_remove:
            os.remove(blastdb_path)
    
    # Get the name of the phage database
    blastdb_path = f"{blastdb_path}.gz"

    all_phages_name = []
    with gzip.open(blastdb_path, "rt") as f:
        parser = SeqIO.parse(f, "fasta")
        for record in parser:
            # Remove the version number from the phage name
            all_phages_name.append(record.id.split('.')[0])

    return all_phages_name


####################################################################################################


def check_mash_index(mash_index_path: str, install: bool = False) -> None:
    """
    Checks if the mash index is present and if not downloads it and creates the index

    Args:
        mash_index_path (str): Path to the mash index
        install (bool, optional): Whether to install the mash index. Defaults to False.

    Returns:
        None
    """

    # check if mash index is present
    if os.path.exists(mash_index_path):
        print_ok(f"Found {mash_index_path} as expected")
    else:
        if install:
            print_error(
                f"File {mash_index_path} does not exist will create database now  "
            )
            print_error("Will download the database now and create database")

            url = "https://millardlab-taxmyphage.s3.climb.ac.uk/ICTV_MSL39v1.msh"

            create_folder(os.path.dirname(mash_index_path))

            # Download the file from the URL to the output directory
            download(url, mash_index_path)
        else:
            print_error(
                f"File {mash_index_path} does not exist. Please download the database and create the database."
                "Or use the --install flag to download and create the database."
            )
            sys.exit()


####################################################################################################

def check_blastn_dataframe(blastn_df_path: str, install: bool = False) -> None:
    """
    Checks if the precomputed blastn dataframe  is present and if not downloads it

    """


    if os.path.exists(blastn_df_path):
        print_ok(f"Found {blastn_df_path} as expected")
    else:
        if install:
            print_error(
                f"File {blastn_df_path} does not exist will create database now  "
            )
            print_error("Will download the database now and create database")

            url = "https://millardlab-taxmyphage.s3.climb.ac.uk/M.pa"

            create_folder(os.path.dirname(blastn_df_path))

            # Download the file from the URL to the output directory
            download(url, blastn_df_path)
        else:
            print_error(
                f"File {blastn_df_path} does not exist. Please download the dataframe or precomputed will not be used ."
                "Or use the --install flag to download and create the dataframe."
            )
            sys.exit()


####################################################################################################


def check_VMR(VMR_path: str, install: bool = False) -> None:
    """
    Checks if the VMR is present and if not downloads it

    Args:
        VMR_path (str): Path to the VMR
        install (bool, optional): Whether to install the VMR. Defaults to False.

    Returns:
        None
    """

    # check if VMR is present
    if os.path.exists(VMR_path):
        print_ok(f"Found {VMR_path} as expected")
    else:
        if install:
            print_error(f"File {VMR_path} does not exist will try downloading now")
            print_error("Will download the current VMR now")

            url = "https://ictv.global/vmr/current"

            create_folder(os.path.dirname(VMR_path))

            # Download the file from the URL to the output directory
            download(url, VMR_path)
        else:
            print_error(
                f"File {VMR_path} does not exist. Please download the database and create the database."
                "Or use the --install flag to download and create the database."
            )
            sys.exit()


####################################################################################################


def install_db(
    VMR_path: str,
    blastdb_path: str,
    mash_index_path: str,
    blastn_df_path: str,
    output: str,
    makeblastdb: str,
) -> None:
    """
    Install the databases and quit

    Args:
        VMR_path (str): Path to the VMR
        blastdb_path (str): Path to the blastDB
        mash_index_path (str): Path to the mash index
        blastn_df_path (str): Path to the precomputed blastn results
        output (str): Path to the output directory
        makeblastdb (str): Path to the makeblastdb executable

    Returns:
        None
    """

    # Check if the VMR file exists
    check_VMR(VMR_path, install=True)

    # Check if the mash index exists
    check_mash_index(mash_index_path, install=True)

    # Check if the mash index exists
    check_blastn_dataframe(blastn_df_path, install=True)
    
    # Download the VMR and create the mash index
    check_blastDB(blastdb_path, output, makeblastdb, install=True)

    print_ok("All databases installed successfully!\n")

    sys.exit()


####################################################################################################
