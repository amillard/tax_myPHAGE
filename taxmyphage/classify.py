"""
This module provides functionalities for classifying bacteriophages based on their genomes.
It includes functions for running mash, parsing mash output, running BLAST, parsing BLAST output,
and performing the classification.
"""

####################################################################################################
# IMPORTS
####################################################################################################

import io
import os
import subprocess
from typing import Tuple, Dict, List
import pandas as pd
from Bio import SeqIO
from icecream import ic

from taxmyphage.pmv import ClusteringOnGenomicSimilarity
from taxmyphage.plot import heatmap
from taxmyphage.utils import print_error, print_ok, print_res, print_warn, create_folder
from taxmyphage.utils import (
    statement_current_genus_new_sp,
    statement_current_genus_sp,
    summary_statement1,
    summary_statement_inconsitent,
    print_table,
)

####################################################################################################
# FUNCTIONS
####################################################################################################

def run_mash(
    query: str, mash_index_path: str, dist: float, threads: int, mash_exe: str
) -> pd.DataFrame:
    """
    Runs mash dist on the query genome against the mash index

    Args:
        query (str): Path to the query genome
        mash_index_path (str): Path to the mash index
        dist (float): Mash distance
        threads (int): Number of threads to use
        mash_exe (str): Path to the mash executable

    Returns:
        pd.DataFrame: Dataframe of the mash results
    """

    # run mash to get top hit and read into a pandas dataframe
    cmd = f"{mash_exe} dist -d {dist} -p {threads} {mash_index_path} {query}"
    #ic(cmd)
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
    mash_exe: str,
    blastdbcmd_exe: str,
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
        mash_exe (str): Path to the mash executable
        blastdbcmd_exe (str): Path to the blastdbcmd executable

    Returns:
        mash_df (pd.DataFrame): Dataframe of the mash results
        accession_genus_dict (Dict[str, str]): Dictionary of accessions linking to genus
    """

    # create a dictionary of Accessions linking to Genus
    accession_genus_dict = taxa_df.set_index("Genbank")["Genus"].to_dict()

    mash_df = run_mash(query, mash_index_path, dist, threads, mash_exe)

    number_hits = mash_df.shape[0]

    # get the number of genomes wih mash distance < args.dist
    if number_hits < 1:
        print_error(
                "Error: No hits were found with the default settings.\n"
                "If this is a phage sequence, it likely represents a new species and genus.\n"
                "However tax_my_phage is unable to classify it at this point as it can only classify at the Genus/Species level.\n"
        )

        with open(taxa_csv_output_path, "w", encoding="utf-8") as wt:
            wt.write("No hits were found with the default settings\n")

        return pd.DataFrame(), accession_genus_dict
    else:
        print_res(f"Number of phage genomes detected with mash distance of < {dist} is:{number_hits}")

    # sort dataframe by distance so they are at the top
    mash_df = mash_df.sort_values(by="distance", ascending=True)

    # write the mash results to a file
    mash_df.to_csv(os.path.join(results_path, "mash.txt"), index=False)

    # print the mash results
    minimum_value = mash_df["distance"].min()
    maximum_value = mash_df["distance"].max()

    print_ok("\nThe mash distances obtained for this query phage:\n")
    print_table(
        {
            "Minimum": [minimum_value],
            "Maximum": [maximum_value],
        },
        print_func=print_ok,
    )

    # set the maximum number of hits to take forward. Max is 10 or the max number in the table if <10
    filter_hits = ""
    if number_hits < 10:
        filter_hits = number_hits
    else:
        filter_hits = 10

    # copy top 10 hits to a new dataframe
    top_10 = mash_df.iloc[:filter_hits].copy()

    # ic(mash_df.head(10))
    # ic(top_10)

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

    ic(top_10[["acc", "genus", "Genus", "Virus name(s)"]].head())

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
        get_genomes_cmd = f"{blastdbcmd_exe} -db {blastdb_path} -entry {','.join(keys)} -out {known_taxa_path}"
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

        get_genomes_cmd = f"{blastdbcmd_exe} -db {blastdb_path} -entry {','.join(list_of_genus_accessions)} -out {known_taxa_path}"
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

def classification_similarity(
    known_taxa_path: str,
    query: str,
    # precomputed code addition
    db_dir:str,
    dont_use_precomputed: bool,
    taxa_df: pd.DataFrame,
    taxa_csv_output_path: str,
    results_path: str,
    threads: int,
    accession_genus_dict: Dict[str, str],
    Figure: bool,
    verbose: bool,
    blastn_exe: str,
    makeblastdb_exe: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    """
    Classifies the query genome using similarity

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
        blastn_exe (str): Path to the blastn executable
        makeblastdb_exe (str): Path to the makeblastdb executable

    Returns:
        merged_df (pd.DataFrame): Dataframe of the merged results of similarity and ICTV dataframe without the query
        query_merged_df pd.DataFrame): Dataframe of the merged results of similarity and ICTV dataframe with the query
        closest_genome (str): Closest genome
    """

    # store files for similarity run- or equivalent
    similarity_folder = os.path.join(results_path, "pmv")

    create_folder(similarity_folder)

    similarity_in_path = os.path.join(similarity_folder, "pmv_in.fa")

    heatmap_file = os.path.join(results_path, "heatmap")
    top_right_matrix = os.path.join(results_path, "top_right_matrix.tsv")
    similarities_file = os.path.join(results_path, "similarities.tsv")

    # Merge the query and known taxa into a single file
    with open(similarity_in_path, "w", encoding="utf-8") as merged_file:
        list_genomes = [known_taxa_path, query]
        for file in list_genomes:
            SeqIO.write(SeqIO.parse(file, "fasta"), merged_file, "fasta")

    # run ClusteringOnGenomicSimilarity
    pmv = ClusteringOnGenomicSimilarity(
        file=similarity_in_path,
        reference=similarity_in_path,
        # precomputed code addition
        db_dir = db_dir,
        dont_use_precomputed = dont_use_precomputed,
        nthreads=threads,
        verbose=verbose,
        blastn_exe=blastn_exe,
        makeblastdb_exe=makeblastdb_exe,
    )
    df1, pmv_outfile = pmv.run()

    # ic(df1)
    # ic(pmv_outfile)
    # ic(pmv.dfM)

    # heatmap and distances
    if Figure:
        print_ok("\nWill calculate and save heatmaps now\n")
        heatmap(pmv.dfM, heatmap_file, top_right_matrix, accession_genus_dict)
    else:
        print_error("\nSkipping calculating heatmaps and saving them\n ")

    pmv.save_similarities(similarities_file)

    # merge the ICTV dataframe with the results of similarity
    # fill in missing with Not Defined yet
    merged_df = pd.merge(
        df1, taxa_df, left_on="genome", right_on="Genbank", how="left"
    ).fillna("Not Defined Yet")

    ic(merged_df.head())

    # write dataframe to file
    merged_df.to_csv(taxa_csv_output_path, sep="\t", index=False)

    # create a copy of this dataframe for later use
    query_merged_df = merged_df.copy()

    merged_df = merged_df[~merged_df["genome"].str.contains("query_")].reset_index(
        drop=True
    )

    # get the closest species
    closest_genome = pmv.get_query_closest()

    return merged_df, query_merged_df, closest_genome


####################################################################################################

def new_genus(
    query_genus_cluster_number: int,
    dict_genus_cluster_2_genus_name: Dict[int, str],
    summary_output_path: str,
    message: str,
) -> Dict[str, str]:
    """
    Classifies the query genome as a new genus

    Args:
        query_genus_cluster_number (int): Query genus cluster number
        dict_genus_cluster_2_genus_name (Dict[int, str]): Dictionary of genus cluster to genus name
        summary_output_path (str): Path to the summary output file
        prefix (str): Prefix to add to the output file
        message (str): Informative message giving the confidence of the annotation

    Returns:
        Dict[str, str]: Dictionary of taxonomic information
    """

    print_warn(
        f"Cluster Number: {query_genus_cluster_number} is not in the dictionary of " 
        f"known Genera: {dict_genus_cluster_2_genus_name}"
    )

    if summary_output_path:
        with open(summary_output_path, "w", encoding="utf-8") as file:
            file.write("--------------- TAXMYPHAGE RESULTS ---------------\n\n")
            file.write(
                "Try running again with if you larger distance if you want a Figure. "
                "The query is both a new genus and species\n\n"
                "Genus: New genus\n"
                "Species: New species\n\n"
                "-------------- INFORMATION MESSAGES --------------\n\n"
                f"INFO: {message}\n"
            )

    return {
        "Realm": "Unknown",
        "Kingdom": "Unknown",
        "Phylum": "Unknown",
        "Class": "Unknown",
        "Order": "Unknown",
        "Family": "Unknown",
        "Subfamily": "Unknown",
        "Genus": "New_genus",
        "Species": "New_species",
        "Message": message,
    }


####################################################################################################


def current_genus_current_species(
    query_species_cluster_number: int,
    dict_species_cluster_2_species_name: Dict[int, str],
    summary_output_path: str,
    dict_genus_cluster_2_genus_name: Dict[int, str],
    query_genus_cluster_number: int,
    merged_df: pd.DataFrame,
    mash_df: pd.DataFrame,
    message: str,
) -> Dict[str, str]:
    """
    Classifies the query genome as a current genus and current species

    Args:
        query_species_cluster_number (int): Query species cluster number
        dict_species_cluster_2_species_name (Dict[int, str]): Dictionary of species cluster to species name
        summary_output_path (str): Path to the summary output file
        dict_genus_cluster_2_genus_name (Dict[int, str]): Dictionary of genus cluster to genus name
        query_genus_cluster_number (int): Query genus cluster number
        merged_df (pd.DataFrame): Dataframe of the merged results of similarity and ICTV dataframe without the query
        mash_df (pd.DataFrame): Dataframe of the mash results
        message (str): Informative message giving the confidence of the annotation

    Returns:
        Dict[str, str]: Dictionary of taxonomic information
    """

    print(
        "Phage is within a current genus and same as a current species\n"
        "....working out which one now .....\n"
    )
    predicted_genus = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
    predicted_species = dict_species_cluster_2_species_name[
        query_species_cluster_number
    ]
    print(
        f"QUERY is in the genus: {predicted_genus} and is species: {predicted_species}"
    )
    # identify the row in the pandas data frame that is the same species
    matching_species_row = merged_df[merged_df["Species"] == predicted_species]
    ic(matching_species_row)

    list_of_S_data = matching_species_row.iloc[0].to_dict()

    list_of_S_data["Message"] = message

    ic(list_of_S_data)

    print_res(
        "\nQuery sequence is:\n"
        f'Class: {list_of_S_data["Class"]}\n'
        f'Family: {list_of_S_data["Family"]}\n'
        f'Subfamily: {list_of_S_data["Subfamily"]}\n'
        f'Genus: {list_of_S_data["Genus"]}\n'
        f'Species: {list_of_S_data["Species"]}\n'
    )

    if summary_output_path:
        with open(summary_output_path, "w", encoding="utf-8") as file:
            file.write("--------------- TAXMYPHAGE RESULTS ---------------\n\n")

            file.write(
                f'{statement_current_genus_sp}'
                f'Class: {list_of_S_data["Class"]}\nFamily: {list_of_S_data["Family"]}\n'
                f'Subfamily: {list_of_S_data["Subfamily"]}\nGenus: {list_of_S_data["Genus"]}\n'
                f'Species: {list_of_S_data["Species"]}\n\n'
                "-------------- INFORMATION MESSAGES --------------\n\n"
                f'INFO: {message}\n'
                f'{summary_statement1}'
            )

            file.write("------------------ MASH RESULTS ------------------\n\n")

        mash_df.to_csv(summary_output_path, mode="a", header=True, index=False, sep="\t")

    return list_of_S_data


####################################################################################################


def current_genus_new_species(
    summary_output_path: str,
    dict_genus_cluster_2_genus_name: Dict[int, str],
    query_genus_cluster_number: int,
    merged_df: pd.DataFrame,
    mash_df: pd.DataFrame,
    message: str,
) -> Dict[str, str]:
    """
    Classifies the query genome as a current genus and new species

    Args:
        summary_output_path (str): Path to the summary output file
        dict_genus_cluster_2_genus_name (Dict[int, str]): Dictionary of genus cluster to genus name
        query_genus_cluster_number (int): Query genus cluster number
        merged_df (pd.DataFrame): Dataframe of the merged results of similarity and ICTV dataframe without the query
        mash_df (pd.DataFrame): Dataframe of the mash results
        message (str): Informative message giving the confidence of the annotation

    Returns:
        Dict[str, str]: Dictionary of taxonomic information
    """

    print(
        "Phage is within a current genus, BUT is representative of a new species\n"
        "....working out which one now .....\n"
    )

    predicted_genus = dict_genus_cluster_2_genus_name[query_genus_cluster_number]

    matching_genus_rows = merged_df[merged_df["Genus"] == predicted_genus]

    dict_exemplar_genus = matching_genus_rows.iloc[0].to_dict()
    genus_value = dict_exemplar_genus["Genus"]
    dict_exemplar_genus["Species"] = genus_value + " new_name"

    dict_exemplar_genus["Message"] = message

    ic(matching_genus_rows)
    ic(genus_value)

    print_res(
        "\nQuery sequence is:\n"
        f'Class: {dict_exemplar_genus["Class"]}\n'
        f'Family: {dict_exemplar_genus["Family"]}\n'
        f'Subfamily: {dict_exemplar_genus["Subfamily"]}\n'
        f'Genus: {dict_exemplar_genus["Genus"]}\n'
        f'Species: {dict_exemplar_genus["Species"]}\n'
    )

    if summary_output_path:
        with open(summary_output_path, "w", encoding="utf-8") as file:
            file.write("--------------- TAXMYPHAGE RESULTS ---------------\n\n")
            
            file.write(
                f'{statement_current_genus_new_sp}'
                f'Class: {dict_exemplar_genus["Class"]}\nFamily: {dict_exemplar_genus["Family"]}\n'
                f'Subfamily: {dict_exemplar_genus["Subfamily"]}\nGenus: {dict_exemplar_genus["Genus"]}\n'
                f'Species: {dict_exemplar_genus["Species"]}\n\n'
                '-------------- INFORMATION MESSAGES --------------\n\n'
                f'INFO: {message}\n'
                f'{summary_statement1}'
            )

            file.write("------------------ MASH RESULTS ------------------\n\n")

        mash_df.to_csv(summary_output_path, mode="a", header=True, index=False, sep="\t")

    return dict_exemplar_genus


####################################################################################################
def new_genus_new_species(
    summary_output_path: str, 
    mash_df: pd.DataFrame,
    message: str,
) -> Dict[str, str]:
    """
    Classifies the query genome as a new genus and new species

    Args:
        summary_output_path (str): Path to the summary output file
        mash_df (pd.DataFrame): Dataframe of the mash results
        message (str): Informative message giving the confidence of the annotation

    Returns:
        Dict[str, str]: Dictionary of taxonomic information
    """

    print(
        "\nQuery does not fall within a current genus or species as defined by ICTV. "
        "Therefore the query sequence is likely the first representative of both a new species and new genus.\n"
        "Data produced by taxmyphage will help you write a Taxonomy proposal so it can be offically classified.\n"
    )

    print_warn(
        "WARNING:: taxmyphage does not compare against all other known phages, only those that have been classified\n"
    )

    if summary_output_path:
        with open(summary_output_path, "w", encoding="utf-8") as file:
            file.write("--------------- TAXMYPHAGE RESULTS ---------------\n\n")

            file.write(
                "Query sequence can not be classified within a current genus or species, it is in.\n\n"
                "Remember taxmyphage compared against viruses classified by the ICTV. Allowing determine if it represents a new "
                "species or geneus.\nIt does not tell you if it is similar to other phages that have yet to be classified \n"
                "You can do this by comparison with INPHARED database if you wish.\n"
            )
        
            file.write("------------------ MASH RESULTS ------------------\n\n")

        mash_df.to_csv(summary_output_path, mode="a", header=True, index=False, sep="\t")

    return {
        "Realm": "Unknown",
        "Kingdom": "Unknown",
        "Phylum": "Unknown",
        "Class": "Unknown",
        "Order": "Unknown",
        "Family": "Unknown",
        "Subfamily": "Unknown",
        "Genus": "New_genus",
        "Species": "New_species",
        "Message": message,
    }


####################################################################################################

def assess_taxonomic_info(
    query_genus_cluster_number: int,
    query_species_cluster_number: int,
    list_ICTV_genus_clusters: List[int],
    list_ICTV_species_clusters: List[int],
    dict_genus_cluster_2_genus_name: Dict[int, str],
    dict_species_cluster_2_species_name: Dict[int, str],
    summary_output_path: str,
    merged_df: pd.DataFrame,
    mash_df: pd.DataFrame,
    message: str,
) -> Dict[str, str]:
    """
    Classifies the query genome as a new genus and new species

    Args:
        query_genus_cluster_number (int): Query genus cluster number
        query_species_cluster_number (int): Query species cluster number
        list_ICTV_genus_clusters (List[int]): List of ICTV genus clusters
        list_ICTV_species_clusters (List[int]): List of ICTV species clusters
        dict_genus_cluster_2_genus_name (Dict[int, str]): Dictionary of genus cluster to genus name
        dict_species_cluster_2_species_name (Dict[int, str]): Dictionary of species cluster to species name
        summary_output_path (str): Path to the summary output file
        merged_df (pd.DataFrame): Dataframe of the merged results of similarity and ICTV dataframe without the query
        mash_df (pd.DataFrame): Dataframe of the mash results
        message (str): Informative message giving the confidence of the annotation

    Returns:
        Dict[str, str]: Dictionary of taxonomic information
    """

    # GENUS CHECK FIRST- Current genus and current species
    if (
        query_genus_cluster_number in list_ICTV_genus_clusters
        and query_species_cluster_number in list_ICTV_species_clusters
    ):
        # print the information that the query already have a genus and a species
        taxonomic_info = current_genus_current_species(
            query_species_cluster_number=query_species_cluster_number,
            dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
            summary_output_path=summary_output_path,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            query_genus_cluster_number=query_genus_cluster_number,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

        # WRITE CODE FOR GIVING INFO ON SPECIES

    # SAME GENUS but different species
    elif (
        query_genus_cluster_number in list_ICTV_genus_clusters
        and query_species_cluster_number not in list_ICTV_species_clusters
    ):
        # print the information that the query already have a genus but not a species
        taxonomic_info = current_genus_new_species(
            summary_output_path=summary_output_path,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            query_genus_cluster_number=query_genus_cluster_number,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

    # NEW GENUS and NEW SPECIES
    elif (
        query_genus_cluster_number not in list_ICTV_genus_clusters
        and query_species_cluster_number not in list_ICTV_species_clusters
    ):
        taxonomic_info = new_genus_new_species(
            summary_output_path=summary_output_path, 
            mash_df=mash_df,
            message=message,
        )

    return taxonomic_info


####################################################################################################

def classification(
    merged_df: pd.DataFrame,
    query_merged_df: pd.DataFrame,
    results_path: str,
    mash_df: pd.DataFrame,
    prefix: str,
    closest_genome: str,
) -> Dict[str, str]:
    """
    Classifies the query genome

    Args:
        merged_df (pd.DataFrame): Dataframe of the merged results of similarity and ICTV dataframe without the query
        query_merged_df (pd.DataFrame): Dataframe of the merged results of similarity and ICTV dataframe with the query
        results_path (str): Path to the results directory
        mash_df (pd.DataFrame): Dataframe of the mash results
        prefix (str): Prefix to add to the output file
        closest_genome (str): Genome name

    Returns:
        Dict[str, str]: Dictionary of taxonomic information
    """

    ic.disable()

    # path the final results summary file
    summary_results = prefix + "Summary_file.txt"
    summary_output_path = os.path.join(results_path, summary_results)

    # Count the number genera
    # excluding query
    num_unique_similarity_genus_clusters = merged_df["genus_cluster"].nunique()
    num_unique_ICTV_genera = merged_df["Genus"].nunique()

    # including query
    total_num_similarity_genus_clusters = query_merged_df["genus_cluster"].nunique()
    total_num_similarity_species_clusters = query_merged_df["species_cluster"].nunique()

    print(
        "\nTotal number of the clustering on genomic similarity algorithm genus clusters in the input including "
        f"QUERY sequence was: {total_num_similarity_genus_clusters}\n"
        "Total number of the clustering on genomic similarity algorithm species clusters including "
        f"QUERY sequence was: {total_num_similarity_species_clusters}"
    )

    print(
        f"\nNumber of current ICTV defined genera was: {num_unique_ICTV_genera}\n"
        f"Number of the clustering on genomic similarity algorithm predicted genera (excluding query) was: {num_unique_similarity_genus_clusters}"
    )

    if num_unique_ICTV_genera == num_unique_similarity_genus_clusters:
        print(
            "\n\nCurrent ICTV and the clustering on genomic similarity algorithm predictions are consistent for the data that was used to compare against"
        )

    print_ok(
        f"\nNumber of unique the clustering on genomic similarity algorithm clusters at default cutoff of 70% is: {num_unique_similarity_genus_clusters}"
    )

    print_ok(
        f"Number of current ICTV genera associated with the reference genomes is: {num_unique_ICTV_genera}"
    )

    # get information on the query from the dataframe
    # get species and genus cluster number
    query_row = query_merged_df[query_merged_df["genome"].str.contains("query_")]

    query_genus_cluster_number = query_row["genus_cluster"].values[0]
    query_species_cluster_number = query_row["species_cluster"].values[0]

    print(f"\nSpecies cluster number is: {query_species_cluster_number}")
    print(f"Genus cluster number is: {query_genus_cluster_number}")

    # list of similarity genus and species numbers
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

    # change the value for the closest genome as the one for the cluster number
    # this is to avoid the problem of having a random name chosen for the group
    # the query is in (e.g. if in the cluster of the query the genus could be
    # Agtrevirus or Limestonevirus, the value in the dictionary will be the one of
    # of the closest genome to the query [using the similarity matrix])

    cluster_closest_genus = (
        merged_df[merged_df["genome"] == closest_genome][["genus_cluster", "Genus"]]
        .values[0]
        .tolist()
    )
    cluster_closest_species = (
        merged_df[merged_df["genome"] == closest_genome][["species_cluster", "Species"]]
        .values[0]
        .tolist()
    )

    ic(cluster_closest_genus)
    ic(cluster_closest_species)

    dict_genus_cluster_2_genus_name[cluster_closest_genus[0]] = cluster_closest_genus[1]
    dict_species_cluster_2_species_name[
        cluster_closest_species[0]
    ] = cluster_closest_species[1]

    ic(dict_genus_cluster_2_genus_name)
    ic(
        sorted(list(dict_genus_cluster_2_genus_name.keys()))
        == sorted(list_ICTV_genus_clusters)
    )
    ic(sorted(list(dict_genus_cluster_2_genus_name.keys())))
    ic(sorted(list(list_ICTV_genus_clusters)))

    # check query is within a current genus. If not, then new Genus
    if query_genus_cluster_number not in dict_genus_cluster_2_genus_name:
        # print the information that the query is a new genus
        message = "Query is a new genus and species. You could try running again with if you larger distance"

        taxonomic_info = new_genus(
            query_genus_cluster_number,
            dict_genus_cluster_2_genus_name,
            summary_output_path,
            message,
        )

        # no more analysis to do so return
        return taxonomic_info

    # get the predicted genus name
    predicted_genus_name = dict_genus_cluster_2_genus_name[query_genus_cluster_number]
    print(f"\nPredicted genus is: {predicted_genus_name}\n")

    # create a dict of species to species_cluster
    # if number of ICTV genera and predicted VIRIDIC-like genera match:
    if num_unique_ICTV_genera == num_unique_similarity_genus_clusters:
        message = "Current ICTV taxonomy and the clustering on genomic similarity algorithm output appear to be consistent at the genus level"
        
        print(message)

        taxonomic_info = assess_taxonomic_info(
            query_genus_cluster_number=query_genus_cluster_number,
            query_species_cluster_number=query_species_cluster_number,
            list_ICTV_genus_clusters=list_ICTV_genus_clusters,
            list_ICTV_species_clusters=list_ICTV_species_clusters,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
            summary_output_path=summary_output_path,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

    # if number of VIRIDIC-like genera is greater than ICTV genera
    elif num_unique_ICTV_genera != num_unique_similarity_genus_clusters:
        print_error(f"""{summary_statement_inconsitent}\n""")

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        taxonomic_info = assess_taxonomic_info(
            query_genus_cluster_number=query_genus_cluster_number,
            query_species_cluster_number=query_species_cluster_number,
            list_ICTV_genus_clusters=list_ICTV_genus_clusters,
            list_ICTV_species_clusters=list_ICTV_species_clusters,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
            summary_output_path=summary_output_path,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

    return taxonomic_info

####################################################################################################
