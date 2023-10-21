# -*- coding: utf-8 -*-

############################################################################################################
# Imports
############################################################################################################

import os
from argparse import ArgumentParser
from tax_myPHAGE import __version__
from tax_myPHAGE.utils import CheckAction

############################################################################################################
# Functions
############################################################################################################


def cli(args=None):
    """
    Command line interface for tax_myPHAGE
    """

    description = """Takes a phage genome as as fasta file and compares against all phage genomes that are currently classified 
         by the ICTV. It does not compare against ALL phage genomes, just classified genomes. Having found the closet related phages 
         it runs the VIRIDIC--algorithm and parses the output to predict the taxonomy of the phage. It is only able to classify to the Genus and Species level"""

    parser = ArgumentParser(
        description=description, conflict_handler="resolve", prog="tax_myPHAGE"
    )

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        help="Show the conda-prefix-replacement version number and exit.",
        version=f"tax_myPHAGE {__version__}",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
    )

    general_option = parser.add_argument_group(title="General options")

    general_option.add_argument(
        "-i",
        "--input",
        dest="in_fasta",
        type=str,
        help="Path to an input fasta file(s), or directory containing fasta files. The fasta file(s) could"
        " contain multiple phage genomes. They can be compressed (gzip). If a directory is given the expected fasta extentions"
        " are ['fasta', 'fna', 'fsa', 'fa'] but can be gzipped. (Required)",
        required=True,
        nargs="+",
    )
    general_option.add_argument(
        "-o",
        "--output",
        type=str,
        default=os.path.join(os.getcwd(), f"taxmyphage_results"),
        dest="output",
        help="Path to the output directory. (Default is current directory)",
    )
    general_option.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="",
        dest="prefix",
        help="will add the prefix to results and summary files that will store results of MASH and comparision to the VMR Data produced by"
        " ICTV combines both sets of this data into a single csv file. "
        " Use this flag if you want to run multiple times and keep the results files without manual renaming of files. (Default no prefix)",
    )
    general_option.add_argument(
        "-d",
        "--distance",
        type=float,
        default=0.2,
        dest="dist",
        help="Will change the mash distance for the intial seraching for close relatives. We suggesting keeping at 0.2"
        " If this results in the phage not being classified, then increasing to 0.3 might result in an output that shows"
        " the phage is a new genus. We have found increasing above 0.2 does not place the query in any current genus, only"
        " provides the output files to demonstrate it falls outside of current genera. (Default is 0.2)",
    )
    general_option.add_argument(
        "--no-figures",
        dest="Figure",
        action="store_false",
        help="Use this option if you don't want to generate Figures. This will speed up the time it takes to run the script"
        " - but you get no Figures. (By default, Figures are generated)",
    )
    general_option.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=str,
        default=1,
        help="Maximum number of threads that will be used by BLASTn. (Default is 1)",
    )

    database_option = parser.add_argument_group(title="Options related to the database")

    database_option.add_argument(
        "-db",
        "--database",
        dest="ICTV_db",
        type=str,
        help="Path to the database of genomes currently classified by the ICTV (If the database does not exist,"
        "it will be downloaded automatically) and stored in the ~/.taxmyPHAGE directory)",
        default=os.path.abspath(
            os.path.join(
                os.path.expanduser("~"), ".taxmyPHAGE", "Bacteriophage_genomes.fasta"
            )
        ),
    )
    database_option.add_argument(
        "--mash_index",
        dest="mash_index",
        type=str,
        help="Path to the prebuilt MASH index of ICTV genomes. (If the index file does not exist,"
        " it will be downloaded automatically and stored in the ~/.taxmyPHAGE directory)",
        default=os.path.abspath(
            os.path.join(os.path.expanduser("~"), ".taxmyPHAGE", "ICTV_2023.msh")
        ),
    )
    database_option.add_argument(
        "--VMR",
        dest="VMR_file",
        type=str,
        help="Path to an input fasta file. (If the VMR file does not exist, it will be downloaded automatically)"
        " and stored in the ~/.taxmyPHAGE directory)",
        default=os.path.abspath(
            os.path.join(os.path.expanduser("~"), ".taxmyPHAGE", "VMR.xlsx")
        ),
    )

    executable_option = parser.add_argument_group(
        title="Executable options, if not in PATH"
    )

    executable_option.add_argument(
        "--blastcmd",
        dest="blastcmd",
        default="blastcmd",
        type=str,
        action=CheckAction,
        help="Path to the blastn executable (default: blastcmd)",
    )
    executable_option.add_argument(
        "--blastn",
        dest="blastn",
        default="blastn",
        type=str,
        action=CheckAction,
        help="Path to the blastn executable (default: blastn)",
    )
    executable_option.add_argument(
        "--makeblastdb",
        dest="makeblastdb",
        default="makeblastdb",
        type=str,
        action=CheckAction,
        help="Path to the blastn executable (default: makeblastdb)",
    )
    executable_option.add_argument(
        "--mashexe",
        dest="mashexe",
        default="mash",
        type=str,
        action=CheckAction,
        help="Path to the MASH executable (default: mash)",
    )

    args, nargs = parser.parse_known_args(args)

    return args, nargs
    # No return value means no error.
    # Return a value of 1 or higher to signify an error.
    # See https://docs.python.org/3/library/sys.html#sys.exit


############################################################################################################

if __name__ == "__main__":
    import sys

    cli(sys.argv[1:])

############################################################################################################
