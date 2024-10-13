"""
This module provides functionalities for handling fasta files in the taxmyphage package.
It includes functions to read and write fasta files, and to count the number of genomes in a fasta file.
"""

####################################################################################################
# Imports
####################################################################################################

import os
import re
import glob
import gzip
from typing import TextIO, List
from Bio import SeqIO
import pandas as pd
from icecream import ic

####################################################################################################
# Functions
####################################################################################################


def read_write_fasta(input_file: str, f: TextIO) -> int:
    """Reads a fasta file and writes it to a new file

    Args:
        input_file (str): Path to the input fasta file
        f (TextIO): File object to write to

    Returns:
        int: Number of genomes in the input file
    """

    handle = (
        gzip.open(input_file, "rt")
        if input_file.endswith(".gz")
        else open(input_file, "rt")
    )

    num = 0

    parser = SeqIO.parse(handle, "fasta")
    for record in parser:
        record.name = record.description = ""
        SeqIO.write(record, f, "fasta")
        num += 1
    handle.close()

    return num


####################################################################################################


def create_files_and_result_paths(
    fasta_files: List[str],
    tmp_fasta: str,
    suffixes: List[str] = ["fasta", "fna", "fsa", "fa"],
) -> int:
    """Creates a multifasta file to parse line by line

    Args:
        fasta_files (List[str]): List of paths to the input fasta files
        tmp_fasta (str): Path to the output fasta file
        suffixes (List[str], optional): List of fasta suffixes. Defaults to ["fasta", "fna", "fsa", "fa"].

    Returns:
        int: Number of genomes in the input file
    """

    fasta_exts = re.compile("|".join([f"\.{suffix}(\.gz)?$" for suffix in suffixes]))
    num_genomes = 0

    with open(tmp_fasta, "w") as f:
        for file in fasta_files:
            if os.path.isdir(file):
                _files = glob.glob(f"{file}/*")
                _files = [x for x in _files if fasta_exts.search(x)]

                for _file in _files:
                    num = read_write_fasta(_file, f)
                    num_genomes += num

            elif os.path.isfile(file):
                num = read_write_fasta(file, f)
                num_genomes += num

    return num_genomes


####################################################################################################


def read_VMR(VMR_path: str) -> pd.DataFrame:
    """Reads the VMR data

    Args:
        VMR_path (str): Path to the VMR data

    Returns:
        pd.DataFrame: DataFrame containing the VMR data
    """
    #update to VMR means sheet posistion changes- change below to 1 or find automatically
    #taxa_df = pd.read_excel(VMR_path, sheet_name=0)
    
    excel_file = pd.ExcelFile(VMR_path)
    sheet_names = excel_file.sheet_names
    # Find the first sheet number that contains 'VMR'
    vmr_sheet_index = next((index for index, sheet in enumerate(sheet_names) if 'VMR' in sheet), None)

    taxa_df = pd.read_excel(VMR_path, sheet_name=vmr_sheet_index)

    # Print the DataFrame and rename a column
    # ic(taxa_df.head())

    # Rename the columns depending on the version of the VMR
    taxa_df = taxa_df.rename(
        columns={"Virus GENBANK accession": "Genbank", "Genome_id": "Genbank"}
    )
    taxa_df["Genbank"] = taxa_df["Genbank"].fillna("")

    return taxa_df


####################################################################################################
