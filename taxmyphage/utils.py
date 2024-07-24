# Desc: Utility functions for tax_myPHAGE

####################################################################################################
# Imports
####################################################################################################

import sys
import os
import shutil
from typing import Dict, Callable

####################################################################################################
# Functions
####################################################################################################


def print_error(txt):
    """
    Print the text in red
    """
    print(f"\033[31m{txt}\033[0m")


def print_warn(txt):
    """
    Print the text in blue (cyan)
    """
    print(f"\033[94m{txt}\033[0m")


def print_ok(txt):
    """
    Print the text in blue
    """
    print(f"\033[34m{txt}\033[0m")


def print_res(txt):
    """
    Print the text in yellow
    """
    print(f"\033[33m{txt}\033[0m")


####################################################################################################


def _validator(cast_func, raw):
    """
    Validate the input
    :param cast_func: function to cast the input
    :type: function
    :param raw: input
    :type: string
    :return: value
    """
    try:
        value = cast_func(raw)
    except ValueError as err:
        print_error(f"Invalid value: {err}")
        sys.exit()
    return value


####################################################################################################


def CheckSoftware(values: str) -> str:
    """
    Check if the software is in the PATH
    :param values: software name
    :type: string
    :return: software name
    """

    def exe(value):
        exe = shutil.which(value)
        if exe:
            return value
        else:
            raise ValueError(f"'{value}' No executable found")

    return _validator(exe, values)


####################################################################################################


def create_folder(mypath):
    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


####################################################################################################

# Statments to output to summary file
summary_statement1 = (
    "\nThe data from the initial mash searching is below as tsv format\n\n"
    "/!\ Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing you determine if "
    "it represents a new species or genus.\nIt does not tell you if it is similar to other phages "
    "that have yet to be classified.\n"
    "You can do this by comparison with INPHARED database if you wish https://github.com/RyanCook94/inphared "
    "or BLAST etc.\n\n"
)


statement_current_genus_new_sp = (
    "Query sequence can be classified within a current genus and represents a new species, it is in:\n\n"
)

statement_current_genus_sp = (
    "Query sequence can be classified within a current genus and species, it is in:\n\n"
)


summary_statement_inconsitent = (
    "The number of expected genera based on current ICTV classification is different from the predicted "
    "number of genus clusters as predicted by clustering on genomic similarity algorithm.\nThis does not mean the current ICTV "
    "classification is wrong (it might be) or that clustering on genomic similarity algorithm is wrong.\nIt could be an edge "
    "case that automated process cannot distinguish. It will require more manual curation to look "
    "at the output files." 
)

####################################################################################################

def print_table(data:Dict[str, float], print_func:Callable[[str], None]=print) -> None:
    """
    Print a formatted table with the given data
    
    Args:
        data (Dict[str, float]): The data to print in the table
        print_func (function, optional): The function to use to print the table. Defaults to print.

    Returns:
        None
    """
    # Get the maximum length of each column
    column_lengths = [max(len(str(value)) for value in column_data) for column_data in data.values()]

    # Print the table header
    header = "+".join(f"+{'-' * (column_length + 2)}" for column_length in column_lengths)
    print_func(f"{header}+")
    print_func(f"|{'|'.join(f' {column_name:^{column_length}} |' for column_name, column_length in zip(data.keys(), column_lengths))}")
    print_func(f"{header}+")

    # Print the table rows
    for row in zip(*data.values()):
        row_str = "|".join(f" {value:^{column_length}} |" if len(str(value)) == column_length else f" {value:^{column_length}} " for value, column_length in zip(row, column_lengths))
        print_func(f"|{row_str}")
        print_func(f"{header}+")

    print()

    return

####################################################################################################
