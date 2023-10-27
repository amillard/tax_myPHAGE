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
from taxmyphage.download_check import check_blastDB, check_mash_index, check_VMR
from taxmyphage.utils import create_folder, print_ok, CheckSoftware
from taxmyphage.handle_files import create_files_and_result_paths, read_VMR
from taxmyphage.classify import (
    classification_mash,
    classification_viridic,
    classification,
)

############################################################################################################

def main():
    # Set up the arguments
    args, nargs = cli.cli()

    for exe in [args.blastdbcmd, args.blastn, args.makeblastdb, args.mash]:
        CheckSoftware(exe)
        
    verbose = args.verbose
    # Defined and set some parameters
    threads = args.threads

    mash_dist = args.dist

    create_folder(args.output)

    # turn on ICECREAM reporting
    if not verbose:
        ic.disable()

    # this is the location of where the script and the databases are stored
    VMR_path = os.path.join(args.db_folder, "VMR.xlsx")
    blastdb_path = os.path.join(args.db_folder, "Bacteriophage_genomes.fasta")
    mash_index_path = os.path.join(args.db_folder, "ICTV_2023.msh")

    print("Looking for database files...\n")

    # Check if the VMR file exists
    check_VMR(VMR_path, install=args.install)

    # read in the VMR data
    taxa_df = read_VMR(VMR_path=VMR_path)

    # Check if the mash index exists
    check_mash_index(mash_index_path, install=args.install)

    # Check if the blastDB file exists
    check_blastDB(blastdb_path, output=args.output, makeblastdb_exe=args.makeblastdb, install=args.install)

    tmp_fasta = os.path.join(args.output, "tmp.fasta")
    # Create a multifasta file to parse line by line

    # create a dictionary to store the taxonomy of each genome
    dict_taxonomy = {}

    num_genomes = create_files_and_result_paths(args.in_fasta, tmp_fasta)

    parser = SeqIO.parse(tmp_fasta, "fasta")

    for genome in tqdm(parser, desc="Classifying", total=num_genomes):

        genome_id = genome.id

        # replace any characters that are not allowed in a file name
        for char in [" ", "/", "|", ":", "(", ")", "[", "]", "{", "}", "<", ">", "#", "%", "&", "+", "$", "="]:
            if char in genome_id:
                genome_id = genome_id.replace(char, "_")

        results_path = os.path.join(args.output, genome_id)
        print_ok(f"\nClassifying {genome.id} in result folder {results_path}...")

        timer_start = time.time()

        print("\nStarting tax_my_phage analysis...\n")

        # create the results folder
        create_folder(results_path)

        # create the path to the query fasta file
        query = os.path.join(results_path, "query.fasta")

        # create a fasta file with just the query genome and add query_ to the id
        with open(query, "w") as output_fid:
            genome.name = genome.description = ""
            genome.id = f"query_{genome_id}"
            SeqIO.write(genome, output_fid, "fasta")

        # path to the combined df containing mash and VMR data
        out_csv_of_taxonomy = args.prefix + "Output_of_taxonomy.csv"
        taxa_csv_output_path = os.path.join(results_path, out_csv_of_taxonomy)

        # fasta file to store known taxa
        known_taxa_path = os.path.join(results_path, "known_taxa.fa")

        mash_df, accession_genus_dict = classification_mash(
            known_taxa_path=known_taxa_path,
            results_path=results_path,
            dist=mash_dist,
            query=query,
            mash_index_path=mash_index_path,
            blastdb_path=blastdb_path,
            taxa_df=taxa_df,
            taxa_csv_output_path=taxa_csv_output_path,
            threads=threads,
            mash_exe=args.mash,
            blastdbcmd_exe=args.blastdbcmd,
        )

        merged_df, copy_merged_df = classification_viridic(
            known_taxa_path=known_taxa_path,
            query=query,
            taxa_df=taxa_df,
            taxa_csv_output_path=taxa_csv_output_path,
            results_path=results_path,
            threads=threads,
            accession_genus_dict=accession_genus_dict,
            Figure=args.Figure,
            verbose=verbose,
            blastn_exe=args.blastn,
            makeblastdb_exe=args.makeblastdb,
        )

        genome_taxo = classification(
            merged_df=merged_df,
            copy_merged_df=copy_merged_df,
            results_path=results_path,
            mash_df=mash_df,
            prefix=args.prefix,
        )

        dict_taxonomy[genome.id] = genome_taxo

        run_time = str(timedelta(seconds=time.time() - timer_start))
        print(f"Run time for {genome.id}: {run_time}\n", file=sys.stderr)
        print("-" * 80, file=sys.stderr)


    # write the taxonomy to a csv file
    taxonomy_tsv = os.path.join(args.output, "Summary_taxonomy.tsv")
    
    with open(taxonomy_tsv, "w") as output_fid:
        output_fid.write("Genome\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tSubfamily\tGenus\tSpecies\tFull_taxonomy\n")
        for key, value in dict_taxonomy.items():
            string_taxo = ""
            full_string = ""
            for taxo in ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species"]:
                # Change the taxonomy to Not Defined Yet if it is empty or nan
                taxonomy = value[taxo] if value[taxo] != "" and value[taxo] == value[taxo] else "Not Defined Yet"
                string_taxo += f"{taxonomy}\t"

                # Change the taxonomy to empty string if it is empty or nan
                taxonomy_full = value[taxo] if value[taxo] != "" and value[taxo] == value[taxo] else ""
                prefix = "sf" if taxo == "Subfamily" else taxo[0].lower()
                full_string += f"{prefix}__{taxonomy_full};"

            # remove the last tab
            string_taxo = string_taxo.rstrip("\t")
            full_string = full_string[:-1]

            # Remove query_ from the genome id
            query = key.replace("query_", "")

            output_fid.write(f"{query}\t{string_taxo}\t{full_string}\n")
            
    # clean up
    os.remove(tmp_fasta)

############################################################################################################

if __name__ == "__main__":
    main()

############################################################################################################
