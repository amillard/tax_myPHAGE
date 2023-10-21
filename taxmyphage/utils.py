# Desc: Utility functions for taxmyphage

############################################################################################################
# Imports
############################################################################################################

import sys, os
import shutil
import argparse
from textwrap import dedent
from icecream import ic

############################################################################################################
# Functions
############################################################################################################


def print_error(txt):
    print(f"\033[31m{txt}\033[0m")


def print_warn(txt):
    print(f"\033[94m{txt}\033[0m")


def print_ok(txt):
    print(f"\033[34m{txt}\033[0m")


def print_res(txt):
    print(f"\033[33m{txt}\033[0m")


############################################################################################################


def _validator(cast_func, raw):
    try:
        value = cast_func(raw)
    except ValueError as err:
        print_error(f"Invalid value: {err}")
        sys.exit()
    return value


############################################################################################################


class CheckAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        def exe(value):
            exe = shutil.which(value)
            if exe:
                return value
            else:
                raise ValueError(f"'{value}' No executable found")

        return _validator(exe, values)


############################################################################################################


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


############################################################################################################

# Statments to output to summary file
summary_statement1 = dedent(
    """
\nThe data from the initial mash searching is below as tsv format\n
Remember taxmyPHAGE compared against viruses classified by the ICTV. Allowing you determine if it represents a new 
species or genus. It does not tell you if it is similar to other phages that have yet to be classified 
You can do this by comparison with INPHARED database if you wish https://github.com/RyanCook94/inphared or BLAST etc \n\n
"""
)

statement_current_genus_new_sp = dedent(
    """
Query sequence can be classified within a current genus and represents a new species, it is in:\n
"""
)

statement_current_genus_sp = dedent(
    """
\nQuery sequence can be classified within a current genus and species, it is in:\n
"""
)

summary_statement_inconsitent = dedent(
    """
The number of expected genera based on current ICTV classification is less than the predicted 
number of genus clusters as predicted by VIRIDIC-algorithm. This does not mean the current ICTV classification
is wrong (it might be)or that VIRIDIC-algorithm is wrong. It could be an edge case that automated process cannot
distinguish. It will require more manual curation to look at the output files.\n 
"""
)

############################################################################################################
