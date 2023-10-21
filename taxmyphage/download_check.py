# Download the VMR and create the mash index

############################################################################################################
# Imports
############################################################################################################

import os
import sys
import subprocess
import wget
import gzip
import shutil

from taxmyphage.utils import print_error, print_ok, create_folder

############################################################################################################
# Functions
############################################################################################################


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


############################################################################################################


def makeblastdb(blastdb_path: str) -> None:
    """
    Creates the blastDB

    Args:
        blastdb_path (str): Path to the blastDB

    Returns:
        None
    """

    makeblastdb_command = f"makeblastdb -in {blastdb_path} -parse_seqids -dbtype nucl"
    try:
        subprocess.run(makeblastdb_command, shell=True, check=True)
        print("makeblastdb command executed successfully!\n")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing makeblastdb: {e}")
        sys.exit()


############################################################################################################


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


############################################################################################################


def check_blastDB(blastdb_path: str, output: str) -> None:
    """
    Checks if the blastDB is present and if not downloads it and creates the database

    Args:
        blastdb_path (str): Path to the blastDB
        output (str): Path to the output directory

    Returns:
        None
    """

    # Variable to check if the uncompressed blastDB needs to be removed to save space
    to_remove = False

    # Check if the blastDB has been created and if not create it
    if os.path.exists(blastdb_path + ".nhr"):
        print_ok(f"Found {blastdb_path}.nhr as expected\n")
    else:
        # check if blastDB is present
        if os.path.exists(blastdb_path):
            print_ok(f"Found {blastdb_path} as expected\n")

            # Unzip the file if it is zipped
            if blastdb_path.endswith(".gz"):
                # Get the basename of the blastdb_path without the .gz extension
                blastdb_path_no_gz = os.path.basename(blastdb_path[:-3])
                blastdb_path_no_gz = os.path.join(output, blastdb_path_no_gz)

                unzip_file(blastdb_path, blastdb_path_no_gz)
                blastdb_path = blastdb_path_no_gz
                to_remove = True
        else:
            print_error(
                f"File {blastdb_path} does not exist will create database now  "
            )
            print_error("Will download the database now and create database")

            url = "https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz"

            create_folder(os.path.dirname(blastdb_path))

            # Download the file from the URL to the output directory
            download(url, f"{blastdb_path}.gz")

            # Gunzip the file
            unzip_file(f"{blastdb_path}.gz", blastdb_path)

            to_remove = True

        makeblastdb(blastdb_path)

        if to_remove:
            os.remove(blastdb_path)


############################################################################################################


def check_mash_index(mash_index_path: str) -> None:
    """
    Checks if the mash index is present and if not downloads it and creates the index

    Args:
        mash_index_path (str): Path to the mash index

    Returns:
        None
    """

    # check if mash index is present
    if os.path.exists(mash_index_path):
        print_ok(f"Found {mash_index_path} as expected\n")
    else:
        print_error(f"File {mash_index_path} does not exist will create database now  ")
        print_error("Will download the database now and create database")

        url = "https://millardlab-inphared.s3.climb.ac.uk/ICTV_2023.msh"

        create_folder(os.path.dirname(mash_index_path))

        # Download the file from the URL to the output directory
        download(url, mash_index_path)


############################################################################################################


def check_VMR(VMR_path: str) -> None:
    """
    Checks if the VMR is present and if not downloads it

    Args:
        VMR_path (str): Path to the VMR

    Returns:
        None
    """

    # check if VMR is present
    if os.path.exists(VMR_path):
        print_ok(f"Found {VMR_path} as expected\n")
    else:
        print_error(f"File {VMR_path} does not exist will try downloading now")
        print_error("Will download the current VMR now")

        url = "https://ictv.global/vmr/current"

        create_folder(os.path.dirname(VMR_path))

        # Download the file from the URL to the output directory
        download(url, VMR_path)


############################################################################################################
