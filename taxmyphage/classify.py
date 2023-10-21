####################################################################################################
# IMPORTS
####################################################################################################

import io
import os
import subprocess
import sys
import time
from datetime import timedelta
from icecream import ic
import pandas as pd
from Bio import SeqIO
from textwrap import dedent
from typing import Tuple, Dict

from taxmyphage.PoorMansViridic import PoorMansViridic
from taxmyphage.plot import heatmap
from taxmyphage.utils import print_error, print_ok, print_res, print_warn
from taxmyphage.utils import (
    statement_current_genus_new_sp,
    statement_current_genus_sp,
    summary_statement1,
    summary_statement_inconsitent,
)

####################################################################################################
# FUNCTIONS
####################################################################################################


def run_mash(
    query: str, mash_index_path: str, dist: float, threads: int
) -> pd.DataFrame:
    """
    Runs mash dist on the query genome against the mash index

    Args:
        query (str): Path to the query genome
        mash_index_path (str): Path to the mash index
        dist (float): Mash distance
        threads (int): Number of threads to use

    Returns:
        pd.DataFrame: Dataframe of the mash results
    """

    # run mash to get top hit and read into a pandas dataframe
    cmd = f"mash dist -d {dist} -p {threads} {mash_index_path} {query}"
    ic(cmd)
    mash_output = subprocess.getoutput(cmd)
    # mash_output = subprocess.check_output(['mash', 'dist', '-d', dist, '-p', threads, mash_index_path, query])

    # list of names for the headers
    mash_df = pd.read_csv(
        io.StringIO(mash_output),
        sep="\t",
        header=None,
        names=["Reference", "Query", "distance", "p-value", "shared-hashes", "ANI"],
    )

    return mash_df


####################################################################################################


def classification_mash(
    known_taxa_path: str,
    results_path: str,
    dist: float,
    query: str,
    mash_index_path: str,
    blastdb_path: str,
    taxa_df: pd.DataFrame,
    taxa_csv_output_path: str,
    threads: int,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Classifies the query genome using mash

    Args:
        known_taxa_path (str): Path to the known taxa
        results_path (str): Path to the results directory
        dist (float): Mash distance
        query (str): Path to the query genome
        mash_index_path (str): Path to the mash index
        blastdb_path (str): Path to the blastDB
        taxa_df (pd.DataFrame): Dataframe of the VMR data
        taxa_csv_output_path (str): Path to the output csv file
        threads (int): Number of threads to use

    Returns:
        mash_df (pd.DataFrame): Dataframe of the mash results
        accession_genus_dict (Dict[str, str]): Dictionary of accessions linking to genus
    """

    # create a dictionary of Accessions linking to Genus
    accession_genus_dict = taxa_df.set_index("Genbank")["Genus"].to_dict()

    mash_df = run_mash(query, mash_index_path, dist, threads)

    number_hits = mash_df.shape[0]

    # get the number of genomes wih mash distance < 0.2
    if number_hits < 1:
        print_error(
            dedent(
                """
                Error: No hits were found with the default settings
                The phage likely represents a new species and genus 
                However tax_my_phage is unable to classify it at this point as it can only classify at the Genus/Species level
                """
            )
        )

        os.system(f"touch {taxa_csv_output_path}")
        sys.exit()
    else:
        print_res(
            dedent(
                f"""
            Number of phage genomes detected with mash distance of < {dist} is:{number_hits}
            """
            )
        )

    # sort dataframe by distance so they are at the top
    mash_df = mash_df.sort_values(by="distance", ascending=True)

    # write the mash results to a file
    mash_df.to_csv(os.path.join(results_path, "mash.txt"), index=False)

    # print the mash results
    minimum_value = mash_df["distance"].min()
    maximum_value = mash_df["distance"].max()

    print_ok(
        dedent(
            f"""\nThe mash distances obtained for this query phage
        is a minimum value of {minimum_value} and maximum value of {maximum_value}\n
        """
        )
    )

    # set the maximum number of hits to take forward. Max is 10 or the max number in the table if <10
    filter_hits = ""
    if number_hits < 10:
        filter_hits = number_hits
    else:
        filter_hits = 10

    # copy top 10 hits to a new dataframe
    top_10 = mash_df.iloc[:filter_hits].copy()

    ic(mash_df.head(10))
    ic(top_10)

    # reindex
    top_10.reset_index(drop=True, inplace=True)

    # get the value at the 10th position
    value_at_10th_position = top_10["distance"].iloc[filter_hits - 1]
    ic(value_at_10th_position)

    # For the top 10 hits, get the genus and accession
    top_10["genus"] = top_10["Reference"].str.split("/").str[1]
    top_10["acc"] = top_10["Reference"].str.split("/").str[-1].str.split(".").str[0]
    top_10 = top_10.merge(taxa_df, left_on="acc", right_on="Genbank")
    top_10["ANI"] = (1 - top_10.distance) * 100

    # returns the unique genera names found in the mash hits - top_10 is not the best name!
    unique_genera_counts = top_10.Genus.value_counts()
    ic(unique_genera_counts.to_dict())
    unique_genera = unique_genera_counts.index.tolist()

    # unique_genera top_10.genus.value_counts().to_dict()
    # print for error checking
    ic(unique_genera)

    # number of genera
    number_of_genera = len(unique_genera)

    print_ok(f"Found {number_of_genera} genera associated with this query genome\n")

    # get all the keys for from a dictionary of accessions and genus names
    keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]

    # print the keys
    ic(keys)

    # Do different things depending how many unique genera were found
    if len(unique_genera) == 1:
        print_ok(
            "Only found 1 genus so will proceed with getting all genomes associated with that genus"
        )

        # get all the keys for from a dictionary of accessions and genus names
        keys = [k for k, v in accession_genus_dict.items() if v == unique_genera[0]]
        number_ok_keys = len(keys)

        print_ok(f"Number of known species in the genus is {number_ok_keys} \n ")

        # create a command string for blastdbcmd
        get_genomes_cmd = f"blastdbcmd -db {blastdb_path} -entry {','.join(keys)} -out {known_taxa_path}"
        res = subprocess.getoutput(get_genomes_cmd)

    elif len(unique_genera) > 1:
        print_ok(
            "Found multiple genera that this query phage might be similar to so will proceed with processing them all"
        )

        list_of_genus_accessions = []
        for i in unique_genera:
            keys = [k for k, v in accession_genus_dict.items() if v == i]
            number_of_keys = len(keys)
            # ic(keys)
            list_of_genus_accessions.extend(keys)
            print_ok(f"Number of known species in the genus {i} is {number_of_keys}")

        ic(list_of_genus_accessions)
        ic(len(list_of_genus_accessions))

        get_genomes_cmd = f"blastdbcmd -db {blastdb_path} -entry {','.join(list_of_genus_accessions)} -out {known_taxa_path}"
        res = subprocess.getoutput(get_genomes_cmd)

    # get smallest mash distance
    min_dist = top_10["distance"].min()

    if min_dist < 0.04:
        print_ok(
            "Phage is likely NOT a new species, will run further analysis now to to confirm this\n"
        )

        top_df = top_10[top_10["distance"] == min_dist]
        ic(top_df)

    elif min_dist >= 0.04 and min_dist < 0.1:
        print_ok(
            "It is not clear if the phage is a new species or not. Will run further analysis now to confirm this...\n"
        )

        top_df = top_10[top_10["distance"] < 0.1]
        ic(top_df)
        ic(top_10.genus.value_counts())

    elif min_dist >= 0.1 and min_dist < dist:
        print_ok("Phage is a new species. Will run further analysis now....\n")

        top_df = top_10[top_10["distance"] >= 0.1]
        ic(top_df)

    return mash_df, accession_genus_dict


####################################################################################################


def classification_viridic(
    known_taxa_path: str,
    query: str,
    taxa_df: pd.DataFrame,
    taxa_csv_output_path: str,
    results_path: str,
    threads: int,
    accession_genus_dict: Dict[str, str],
    Figure: bool,
    verbose: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Classifies the query genome using viridic

    Args:
        known_taxa_path (str): Path to the known taxa
        query (str): Path to the query genome
        taxa_df (pd.DataFrame): Dataframe of the VMR data
        taxa_csv_output_path (str): Path to the output csv file
        results_path (str): Path to the results directory
        threads (int): Number of threads to use
        accession_genus_dict (Dict[str, str]): Dictionary of accessions linking to genus
        Figure (bool): Whether to generate figures
        verbose (bool): Whether to print verbose output

    Returns:
        merged_df (pd.DataFrame): Dataframe of the merged results of VIRIDIC and ICTV dataframe without the query
        copy_merged_df pd.DataFrame): Dataframe of the merged results of VIRIDIC and ICTV dataframe with the query
    """

    # store files for VIRIDIC run- or equivalent
    viridic_in_path = os.path.join(results_path, "viridic_in.fa")

    heatmap_file = os.path.join(results_path, "heatmap")
    top_right_matrix = os.path.join(results_path, "top_right_matrix.tsv")
    similarities_file = os.path.join(results_path, "similarities.tsv")

    # Merge the query and known taxa into a single file
    with open(viridic_in_path, "w") as merged_file:
        list_genomes = [known_taxa_path, query]
        for file in list_genomes:
            SeqIO.write(SeqIO.parse(file, "fasta"), merged_file, "fasta")

    # run VIRIDIC
    PMV = PoorMansViridic(viridic_in_path, nthreads=threads, verbose=verbose)
    df1, pmv_outfile = PMV.run()

    ic(df1)
    ic(pmv_outfile)
    ic(PMV.dfM)

    # heatmap and distances
    if Figure:
        print_ok("\nWill calculate and save heatmaps now\n")
        heatmap(PMV.dfM, heatmap_file, top_right_matrix, accession_genus_dict)
    else:
        print_error("\nSkipping calculating heatmaps and saving them\n ")

    PMV.save_similarities(similarities_file)

    # merge the ICTV dataframe with the results of viridic
    # fill in missing with Not Defined yet
    merged_df = pd.merge(
        df1, taxa_df, left_on="genome", right_on="Genbank", how="left"
    ).fillna("Not Defined Yet")

    ic(merged_df.head())

    # write dataframe to file
    merged_df.to_csv(taxa_csv_output_path, sep="\t", index=False)

    # create a copy of this dataframe for later use
    copy_merged_df = merged_df.copy()

    merged_df = merged_df[~merged_df["genome"].str.contains("query_")].reset_index(
        drop=True
    )

    return merged_df, copy_merged_df


####################################################################################################


def new_genus(
    query_genus_cluster_number: int,
    dict_genus_cluster_2_genus_name: Dict[int, str],
    summary_output_path: str,
    prefix: str,
) -> None:
    """
    Classifies the query genome as a new genus

    Args:
        query_genus_cluster_number (int): Query genus cluster number
        dict_genus_cluster_2_genus_name (Dict[int, str]): Dictionary of genus cluster to genus name
        summary_output_path (str): Path to the summary output file
        prefix (str): Prefix to add to the output file

    Returns:
        None
    """

    print_warn(
        dedent(
            f"""Cluster Number: {query_genus_cluster_number} is not in the dictionary of 
                known Genera: {dict_genus_cluster_2_genus_name}"""
        )
    )

    print_res(
        dedent(
            """
            Phage is NOT within a current genus or species and therefore a both 
            a new Genus and species.\n"""
        )
    )

    with open(summary_output_path, "w") as file:
        file.write(
            f"""Try running again with if you larger distance if you want a Figure.
        The query is both a new genus and species\n
        {prefix}\tNew genus\tNew species\n"""
        )

    return


####################################################################################################


def current_genus_current_species(
    query_species_cluster_number: int,
    dict_species_cluster_2_species_name: Dict[int, str],
    summary_output_path: str,
    dict_genus_cluster_2_genus_name: Dict[int, str],
    query_genus_cluster_number: int,
    merged_df: pd.DataFrame,
    mash_df: pd.DataFrame,
) -> None:
    """
    Classifies the query genome as a current genus and current species

    Args:
        query_species_cluster_number (int): Query species cluster number
        dict_species_cluster_2_species_name (Dict[int, str]): Dictionary of species cluster to species name
        summary_output_path (str): Path to the summary output file
        dict_genus_cluster_2_genus_name (Dict[int, str]): Dictionary of genus cluster to genus name
        query_genus_cluster_number (int): Query genus cluster number
        merged_df (pd.DataFrame): Dataframe of the merged results of VIRIDIC and ICTV dataframe without the query
        mash_df (pd.DataFrame): Dataframe of the mash results

    Returns:
        None
    """

    print(
        dedent(
            """\nPhage is within a current genus and same as a current species 
                ....working out which one now .....\n"""
        )
    )
    predicted_genus = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
    predicted_species = dict_species_cluster_2_species_name[
        query_species_cluster_number
    ]
    print(
        f"""QUERY is in the genus: {predicted_genus} and is species: {predicted_species}"""
    )
    # identify the row in the pandas data frame that is the same species
    matching_species_row = merged_df[merged_df["Species"] == predicted_species]
    ic(matching_species_row)

    list_of_S_data = matching_species_row.iloc[0].to_dict()
    ic(list_of_S_data)

    print_res(
        dedent(
            f"""\nQuery sequence is: 
            Class: {list_of_S_data["Class"]}
            Family: {list_of_S_data["Family"]}
            Subfamily: {list_of_S_data["Subfamily"]}
            Genus: {list_of_S_data["Genus"]}
            Species: {list_of_S_data["Species"]}
            """
        )
    )

    with open(summary_output_path, "w") as file:
        file.write(
            dedent(
                f"""
                {statement_current_genus_sp}\n
                Class: {list_of_S_data["Class"]}\tFamily: {list_of_S_data["Family"]}\t
                Subfamily: {list_of_S_data["Subfamily"]}\tGenus: {list_of_S_data["Genus"]}\t
                Species: {list_of_S_data["Species"]}\n
                {summary_statement1}
                """
            )
        )

    mash_df.to_csv(summary_output_path, mode="a", header=True, index=False, sep="\t")

    return


####################################################################################################


def current_genus_new_species(
    summary_output_path: str,
    query_genus_cluster_number: int,
    merged_df: pd.DataFrame,
    mash_df: pd.DataFrame,
) -> None:
    """
    Classifies the query genome as a current genus and new species

    Args:
        summary_output_path (str): Path to the summary output file
        query_genus_cluster_number (int): Query genus cluster number
        merged_df (pd.DataFrame): Dataframe of the merged results of VIRIDIC and ICTV dataframe without the query
        mash_df (pd.DataFrame): Dataframe of the mash results

    Returns:
        None
    """

    print(
        dedent(
            """\nPhage is within a current genus, BUT is representative of a new species 
            ....working out which one now .....\n"""
        )
    )

    matching_genus_rows = merged_df[
        merged_df["genus_cluster"] == query_genus_cluster_number
    ]

    dict_exemplar_genus = matching_genus_rows.iloc[0].to_dict()
    genus_value = dict_exemplar_genus["Genus"]
    ic(matching_genus_rows)
    ic(genus_value)

    print_res(
        dedent(
            f"""\n
            Query sequence is: 
            Class: {dict_exemplar_genus['Class']}
            Family: {dict_exemplar_genus['Family']} 
            Subfamily: {dict_exemplar_genus['Subfamily']}
            Genus: {dict_exemplar_genus['Genus']}
            Species: {dict_exemplar_genus['Genus']} new_name
            """
        )
    )

    with open(summary_output_path, "a") as file:
        file.write(
            dedent(
                f"""
                {statement_current_genus_new_sp}\n
                Class: {dict_exemplar_genus['Class']}\tFamily: {dict_exemplar_genus['Family']}\t
                Subfamily: {dict_exemplar_genus['Subfamily']}\tGenus: {dict_exemplar_genus['Genus']}\t
                Species: new_specices_name\n
                {summary_statement1}
                """
            )
        )

    mash_df.to_csv(summary_output_path, mode="a", header=True, index=False, sep="\t")


####################################################################################################


def new_genus_new_species(summary_output_path: str, mash_df: pd.DataFrame) -> None:
    """
    Classifies the query genome as a new genus and new species

    Args:
        summary_output_path (str): Path to the summary output file
        mash_df (pd.DataFrame): Dataframe of the mash results

    Returns:
        None
    """

    print(
        dedent(
            """\nQuery does not fall within a current genus or species as defined by ICTV
            Therefore the query sequence is likely the first representative of both a new species and new genus.\n
            Data produced by taxmyphage will help you write a Taxonomy proposal so it can be offically classified.
            \n"""
        )
    )

    print_warn(
        "WARNING:: taxmyphage does not compare against all other known phages, only those that have been classified\n"
    )

    with open(summary_output_path, "w") as file:
        file.write(
            dedent(
                """
                Query sequence can not be classified within a current genus or species, it is in:\n
                Remember taxmyphage compared against viruses classified by the ICTV. Allowing determine if it represents a new 
                species or geneus. It does not tell you if it is similar to other phages that have yet to be classified
                You can do this by comparison with INPHARED database if you wish.\n
                """
            )
        )
    mash_df.to_csv(summary_output_path, mode="a", header=True, index=False, sep="\t")


####################################################################################################


def classification(
    merged_df: pd.DataFrame,
    copy_merged_df: pd.DataFrame,
    results_path: str,
    mash_df: pd.DataFrame,
    prefix: str,
    genome: str,
    timer_start: float,
) -> None:
    """
    Classifies the query genome

    Args:
        merged_df (pd.DataFrame): Dataframe of the merged results of VIRIDIC and ICTV dataframe without the query
        copy_merged_df (pd.DataFrame): Dataframe of the merged results of VIRIDIC and ICTV dataframe with the query
        results_path (str): Path to the results directory
        mash_df (pd.DataFrame): Dataframe of the mash results
        prefix (str): Prefix to add to the output file
        genome (str): Genome name
        timer_start (float): Start time

    Returns:
        None
    """

    # path the final results summary file
    summary_results = prefix + "Summary_file.txt"
    summary_output_path = os.path.join(results_path, summary_results)

    # Count the number genera
    # excluding query
    num_unique_viridic_genus_clusters = merged_df["genus_cluster"].nunique()
    num_unique_ICTV_genera = merged_df["Genus"].nunique()

    # including query
    total_num_viridic_genus_clusters = copy_merged_df["genus_cluster"].nunique()
    total_num_viridic_species_clusters = copy_merged_df["species_cluster"].nunique()

    print(
        dedent(
            f"""\n\nTotal number of VIRIDIC-algorithm genus clusters in the input including QUERY sequence was: {total_num_viridic_genus_clusters}
            Total number of VIRIDIC-algorithm species clusters including QUERY sequence was {total_num_viridic_species_clusters}"""
        )
    )

    print(
        dedent(
            f"""\nNumber of current ICTV defined genera was: {num_unique_ICTV_genera}
            Number of VIRIDIC-algorithm predicted genera (excluding query) was: {num_unique_viridic_genus_clusters}"""
        )
    )

    if num_unique_ICTV_genera == num_unique_viridic_genus_clusters:
        print(
            f"""\n\nCurrent ICTV and VIRIDIC-algorithm predictions are consistent for the data that was used to compare against"""
        )

    print_ok(
        f"\nNumber of unique VIRIDIC-algorithm clusters at default cutoff of 70% is: {num_unique_viridic_genus_clusters}"
    )

    print_ok(
        f"""Number of current ICTV genera associated with the reference genomes is {num_unique_ICTV_genera}"""
    )

    species_genus_dict = merged_df.set_index("species_cluster")["Species"].to_dict()
    ic(species_genus_dict)

    # get information on the query from the dataframe
    # get species and genus cluster number
    query_row = copy_merged_df[copy_merged_df["genome"].str.contains("query_")]

    query_genus_cluster_number = query_row["genus_cluster"].values[0]
    query_species_cluster_number = query_row["species_cluster"].values[0]

    print(
        f"\nCluster number of species is {query_species_cluster_number} and cluster of genus is {query_genus_cluster_number}"
    )

    print(f"Genus cluster number is {query_genus_cluster_number}")

    # list of VIRIDIC genus and species numbers
    list_ICTV_genus_clusters = merged_df["genus_cluster"].unique().tolist()
    list_ICTV_species_clusters = merged_df["species_cluster"].unique().tolist()

    ic(list_ICTV_genus_clusters)
    ic(list_ICTV_species_clusters)

    # create a dictionary linking genus_cluster to genus data
    dict_genus_cluster_2_genus_name = merged_df.set_index("genus_cluster")[
        "Genus"
    ].to_dict()
    dict_species_cluster_2_species_name = merged_df.set_index("species_cluster")[
        "Species"
    ].to_dict()

    ic(dict_genus_cluster_2_genus_name)

    # check query is within a current genus. If not, then new Genus
    if query_genus_cluster_number not in dict_genus_cluster_2_genus_name:
        # print the information that the query is a new genus
        new_genus(
            query_genus_cluster_number,
            dict_genus_cluster_2_genus_name,
            summary_output_path,
            prefix,
        )

        run_time = str(timedelta(seconds=time.time() - timer_start))
        print(f"Run time for {genome.id}: {run_time}\n")
        print("-" * 80)

        # no more analysis to do so return
        return

    # get the predicted genus name
    predicted_genus_name = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
    print(f"\nPredicted genus is: {predicted_genus_name}\n")

    # create a dict of species to species_cluster
    # if number of ICTV genera and predicted VIRIDIC genera match:
    if num_unique_ICTV_genera == num_unique_viridic_genus_clusters:
        print(
            """Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level"""
        )

        # GENUS CHECK FIRST- Current genus and current species
        if (
            query_genus_cluster_number in list_ICTV_genus_clusters
            and query_species_cluster_number in list_ICTV_species_clusters
        ):
            # print the information that the query already have a genus and a species
            current_genus_current_species(
                query_species_cluster_number,
                dict_species_cluster_2_species_name,
                summary_output_path,
                dict_genus_cluster_2_genus_name,
                query_genus_cluster_number,
                merged_df,
                mash_df,
            )

            # WRITE CODE FOR GIVING INFO ON SPECIES

        # SAME GENUS but different species
        elif (
            query_genus_cluster_number in list_ICTV_genus_clusters
            and query_species_cluster_number not in list_ICTV_species_clusters
        ):
            current_genus_new_species(
                summary_output_path, query_genus_cluster_number, merged_df, mash_df
            )

        # NEW GENUS and NEW SPECIES
        elif (
            query_genus_cluster_number not in list_ICTV_genus_clusters
            and query_species_cluster_number not in list_ICTV_species_clusters
        ):
            new_genus_new_species(summary_output_path, mash_df)

    # if number of VIRIDIC genera is greater than ICTV genera
    elif num_unique_ICTV_genera < num_unique_viridic_genus_clusters:
        print_error(f"""{summary_statement_inconsitent}\n""")

        # GENUS CHECK FIRST- Current genus and current species
        if (
            query_genus_cluster_number in list_ICTV_genus_clusters
            and query_species_cluster_number in list_ICTV_species_clusters
        ):
            # print the information that the query already have a genus and a species
            current_genus_current_species(
                query_species_cluster_number,
                dict_species_cluster_2_species_name,
                summary_output_path,
                dict_genus_cluster_2_genus_name,
                query_genus_cluster_number,
                merged_df,
                mash_df,
            )

        # SAME GENUS but different species
        elif (
            query_genus_cluster_number in list_ICTV_genus_clusters
            and query_species_cluster_number not in list_ICTV_species_clusters
        ):
            current_genus_new_species(
                summary_output_path, query_genus_cluster_number, merged_df, mash_df
            )

        # NEW GENUS and NEW SPECIES
        elif (
            query_genus_cluster_number not in list_ICTV_genus_clusters
            and query_species_cluster_number not in list_ICTV_species_clusters
        ):
            new_genus_new_species(summary_output_path, mash_df)

    run_time = str(timedelta(seconds=time.time() - timer_start))
    print(f"Run time for {genome.id}: {run_time}\n", file=sys.stderr)
    print("-" * 80, file=sys.stderr)
