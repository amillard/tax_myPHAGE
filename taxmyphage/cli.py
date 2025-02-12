# -*- coding: utf-8 -*-

"""
This module is used to handle command line interface (CLI) for the taxmyphage package.
It includes classes and functions to parse command line arguments and display help messages.
"""

####################################################################################################
# Imports
####################################################################################################

import os
from argparse import ArgumentParser, ArgumentError
import re
from taxmyphage import __version__

####################################################################################################
# Classes
####################################################################################################


class ArgumentErrorPatch(ArgumentError):
    """ArgumentError that removes subparser metavar from argument choices"""

    def __str__(self):
        msg = super().__str__()
        return (
            re.sub(r" \(choose from .*\)", "", super().__str__())
            if bool(re.match(r"argument \{.*\}", msg))
            else msg
        )


####################################################################################################
# Functions
####################################################################################################


def _metavar(parser, hidden_cmds=set()):
    """Set metavar for subparsers."""
    parser.metavar = (
        "{"
        + ",".join(cmd for cmd in parser._name_parser_map if not cmd in hidden_cmds)
        + "}"
    )


####################################################################################################


def cli(args=None):
    """
    Command line interface for taxmyphage
    """

    description = """Takes a phage genome as as fasta file and compares against all phage genomes that are currently classified
                by the ICTV. It does not compare against ALL phage genomes, just classified genomes. Having found the closet related phages
                it runs the clustering on genomic similarity algorithm and parses the output to predict the taxonomy of the phage. It is only able to classify
                to the Genus and Species level"""

    ArgumentError = ArgumentErrorPatch

    parser = ArgumentParser(
        description=description, conflict_handler="resolve", prog="taxmyphage"
    )

    # Create subparsers
    subparsers = parser.add_subparsers(title="Commands", required=True, dest="command")

    ####################################################################################################
    # Create general subparser that will be given to all subparsers
    ####################################################################################################

    general_subparser = subparsers.add_parser("general", add_help=False)

    general_subparser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show verbose output. (For debugging purposes)",
    )

    general_subparser.add_argument(
        "-V",
        "--version",
        action="version",
        help="Show the version number and exit.",
        version=f"taxmyphage v{__version__}",
    )

    database_option = general_subparser.add_argument_group(title="Databases options")

    database_option.add_argument(
        "-db",
        "--db_folder",
        dest="db_folder",
        metavar="FOLDER_PATH",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "database"),
        help=f"Path to the database directory where the databases are stored. (Default is {os.path.join(os.path.dirname(__file__), 'database')})",
    )

    ####################################################################################################
    # Install subparser
    ####################################################################################################

    install_parser = subparsers.add_parser(
        "install",
        help="Install databases",
        conflict_handler="resolve",
        parents=[general_subparser],
    )

    install_option = install_parser.add_argument_group(title="Install options")

    install_option.add_argument(
        "--makeblastdb",
        dest="makeblastdb",
        default="makeblastdb",
        type=str,
        help="Path to the blastn executable (default: makeblastdb)",
    )

    ####################################################################################################
    # I/O subparser
    ####################################################################################################

    io_parser = subparsers.add_parser("io", add_help=False)

    general_option = io_parser.add_argument_group(title="General options")

    general_option.add_argument(
        "-i",
        "--input",
        dest="in_fasta",
        metavar="[FASTA_FILE ...]",
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
        "-f",
        "--force",
        action="store_true",
        dest="force_overwrite",
        help="Overwrites the genome output directory",
    )
    general_option.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="",
        dest="prefix",
        help="Will add the prefix to results and summary files that will store results of MASH and comparision to the VMR Data produced by"
        " ICTV combines both sets of this data into a single csv file. "
        " Use this flag if you want to run multiple times and keep the results files without manual renaming of files. (Default no prefix)",
    )
    general_option.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        help="Maximum number of threads that will be used by BLASTn. (Default is 1)",
    )
    general_option.add_argument(
        "--no-precomputed",
        action="store_true",
        dest="no_precomputed",
        help="Don't use the precomputed blastn matrix",
    )

    ####################################################################################################
    # MASH options
    ####################################################################################################

    mash_parser_options = subparsers.add_parser("mash_options", add_help=False)

    mash_option = mash_parser_options.add_argument_group(title="MASH options")

    mash_option.add_argument(
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
    mash_option.add_argument(
        "--mash",
        dest="mash",
        default="mash",
        type=str,
        help="Path to the MASH executable (default: mash)",
    )
    mash_option.add_argument(
        "--blastdbcmd",
        dest="blastdbcmd",
        default="blastdbcmd",
        type=str,
        help="Path to the blastn executable (default: blastdbcmd)",
    )

    ####################################################################################################
    # MASH subparser
    ####################################################################################################

    mash_parser = subparsers.add_parser(
        "mash",
        help="Run MASH",
        conflict_handler="resolve",
        parents=[io_parser, mash_parser_options, general_subparser],
    )

    ####################################################################################################
    # Clustering on Genomic Similarity options
    ####################################################################################################

    PMV_parser_options = subparsers.add_parser("genomic_similarity_options", add_help=False)

    comparison_option = PMV_parser_options.add_argument_group(
        title="Comparison options"
    )
    comparison_option.add_argument(
        "--reference",
        dest="reference",
        default="",
        type=str,
        help="Path to the reference database file. Input file will be used as query against it."
        " If not provided, input will be compare against itself."
        " If you use reference no figure is generated. (Default is '')",
    )

    PMV_option = PMV_parser_options.add_argument_group(
        title="Similarity options"
    )

    PMV_option.add_argument(
        "--blastn",
        dest="blastn",
        default="blastn",
        type=str,
        help="Path to the blastn executable (default: blastn)",
    )
    PMV_option.add_argument(
        "--makeblastdb",
        dest="makeblastdb",
        default="makeblastdb",
        type=str,
        help="Path to the blastn executable (default: makeblastdb)",
    )
    PMV_option.add_argument(
        "--no-figures",
        dest="Figure",
        action="store_false",
        help="Use this option if you don't want to generate Figures. This will speed up the time it takes to run the script"
        " - but you get no Figures. (By default, Figures are generated)",
    )

    ####################################################################################################
    # Clustering on Genomic Similarity subparser
    ####################################################################################################

    PMV_parser = subparsers.add_parser(
        "similarity",
        help="Run clustering on genomic similarity",
        conflict_handler="resolve",
        parents=[io_parser, PMV_parser_options, general_subparser],
    )

    ####################################################################################################
    # Run subparser
    ####################################################################################################

    run_parser = subparsers.add_parser(
        "run",
        help="Run taxmyphage",
        conflict_handler="resolve",
        parents=[
            io_parser,
            mash_parser_options,
            PMV_parser_options,
            general_subparser,
        ],
    )

    _metavar(
        subparsers,
        hidden_cmds={"general", "io", "mash_options", "similarity_options", "mash"},
    )

    args, nargs = parser.parse_known_args(args)

    return args, nargs
    # No return value means no error.
    # Return a value of 1 or higher to signify an error.
    # See https://docs.python.org/3/library/sys.html#sys.exit


####################################################################################################
