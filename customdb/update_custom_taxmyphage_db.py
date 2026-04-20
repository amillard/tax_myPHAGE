#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from build_custom_taxmyphage_db import (
    FASTA_EXTENSIONS,
    VMR_COLUMNS,
    build_precomputed_blast_parquet,
    build_summary_tables,
    gzip_copy,
    make_blast_db,
    make_mash_sketch,
    stage_mash_references,
    write_pairwise_similarity_table,
    write_readme,
    write_representatives_fasta,
)

# This script maintains an existing custom taxmyphage database.
#
# The updater does not trust taxmyphage precomputed classification for novel
# incoming genomes. Instead, it recomputes a fresh combined similarity graph
# across:
# - genomes already known to the custom DB
# - one or more new candidate genomes
#
# From that combined graph it reassigns genus/species clusters using the same
# 70/95 thresholds, updates synthetic names, selects representatives, and then
# rebuilds the database assets taxmyphage expects.


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare one or more genomes to an existing custom taxmyphage database, "
            "detect new species/genera, and rebuild the database if needed."
        )
    )
    parser.add_argument(
        "--db-dir",
        type=Path,
        default=Path("custom_taxmyphage_db"),
        help="Existing custom database directory.",
    )
    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="Input FASTA files or directories containing FASTA files.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for BLAST and mash sketch.",
    )
    parser.add_argument(
        "--output-db-dir",
        type=Path,
        default=None,
        help="Output database directory. Default: update --db-dir in place.",
    )
    parser.add_argument(
        "--prefix",
        default="Custom",
        help="Prefix used when assigning new fake genus and species names.",
    )
    parser.add_argument(
        "--keep-workdir",
        action="store_true",
        help="Keep the temporary similarity work directory for inspection.",
    )
    parser.add_argument(
        "--build-mpa",
        action="store_true",
        help="Build M.pa and Mpa_pairwise_similarity.tsv. Default: off for safer no-precomputed querying.",
    )
    return parser.parse_args()


def find_input_fastas(inputs: list[str]) -> list[Path]:
    # Accept either explicit FASTA files or directories of FASTA files.
    fasta_files: list[Path] = []
    for raw_input in inputs:
        input_path = Path(raw_input)
        if input_path.is_dir():
            fasta_files.extend(
                path
                for path in sorted(input_path.iterdir())
                if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
            )
        elif input_path.is_file():
            fasta_files.append(input_path)
        else:
            raise FileNotFoundError(f"Input path not found: {input_path}")
    fasta_files = sorted(dict.fromkeys(path.resolve() for path in fasta_files))
    if not fasta_files:
        raise FileNotFoundError("No FASTA inputs were found")
    return fasta_files


def load_records_from_fastas(fasta_files: list[Path]) -> dict[str, object]:
    # Normalize IDs to filename stems so all downstream tables and FASTA headers
    # use the same stable genome identifier.
    records = {}
    for fasta_path in fasta_files:
        entries = list(SeqIO.parse(str(fasta_path), "fasta"))
        if len(entries) != 1:
            raise ValueError(
                f"{fasta_path} contains {len(entries)} FASTA records. "
                "Expected one genome per FASTA file."
            )
        record = entries[0]
        genome_id = fasta_path.stem
        if genome_id in records:
            raise ValueError(f"Duplicate genome ID detected: {genome_id}")
        record.id = genome_id
        record.name = genome_id
        record.description = genome_id
        records[genome_id] = record
    return records


def load_existing_db_records(db_dir: Path) -> dict[str, object]:
    # all_genomes.fasta is the canonical source for the maintained dataset.
    # For older databases that only stored representatives, fall back to
    # representatives.fasta so the updater can still bootstrap itself.
    all_genomes_fasta = db_dir / "all_genomes.fasta"
    representative_fasta = db_dir / "representatives.fasta"
    source_fasta = all_genomes_fasta if all_genomes_fasta.exists() else representative_fasta
    if not source_fasta.exists():
        raise FileNotFoundError(
            f"Could not find {all_genomes_fasta} or {representative_fasta}"
        )
    records = {}
    for record in SeqIO.parse(str(source_fasta), "fasta"):
        genome_id = record.id.split()[0]
        record.id = genome_id
        record.name = genome_id
        record.description = genome_id
        records[genome_id] = record
    return records


def load_existing_metadata(db_dir: Path, prefix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    # all_genomes_metadata.tsv is the updater's canonical per-genome metadata.
    # If it is missing, reconstruct enough information from the older
    # representative-only files to keep the workflow working.
    all_metadata_path = db_dir / "all_genomes_metadata.tsv"
    raw_vmr_df = read_vmr_rows(db_dir / "VMR.xlsx")
    normalized_vmr_df = raw_vmr_df.rename(columns={"Virus GENBANK accession": "genome"})
    if all_metadata_path.exists():
        metadata_df = pd.read_csv(all_metadata_path, sep="\t")
    else:
        representatives_df = pd.read_csv(db_dir / "representatives.tsv", sep="\t")
        metadata_df = representatives_df.merge(
            normalized_vmr_df[["genome", "Genus", "Species"]],
            on="genome",
            how="left",
        )
        metadata_df["is_representative"] = True
        metadata_df["status"] = "existing_database_member"
        metadata_df["family_name"] = metadata_df["Genus"].apply(
            lambda genus_name: family_from_genus_name(genus_name, prefix)
        )
        metadata_df = metadata_df.rename(
            columns={"Genus": "genus_name", "Species": "species_name"}
        )
    required_columns = {
        "genome",
        "genus_name",
        "species_name",
        "is_representative",
        "status",
        "family_name",
    }
    missing = required_columns - set(metadata_df.columns)
    if missing:
        raise ValueError(f"Existing metadata is missing columns: {sorted(missing)}")
    return metadata_df, raw_vmr_df


def read_vmr_rows(vmr_path: Path) -> pd.DataFrame:
    excel_file = pd.ExcelFile(vmr_path)
    vmr_sheet_index = next(
        (index for index, sheet in enumerate(excel_file.sheet_names) if "VMR" in sheet),
        None,
    )
    if vmr_sheet_index is None:
        raise ValueError(f"Could not find a VMR sheet in {vmr_path}")
    return pd.read_excel(vmr_path, sheet_name=vmr_sheet_index)


def family_from_genus_name(genus_name: str, prefix: str) -> str:
    # Synthetic family names are tied to the custom genus numbering convention.
    match = re.fullmatch(rf"{re.escape(prefix)}_genus_(\d+)", str(genus_name))
    if match:
        return f"{prefix}_family_{int(match.group(1)):03d}"
    return f"{prefix}_family_unassigned"


def choose_canonical_name(names: set[str], prefix: str, label: str) -> str:
    # If an updated clustering merges two older custom groups, keep the
    # lowest-numbered compatible existing name as the canonical one.
    if not names:
        raise ValueError(f"No names supplied for canonical {label} selection")
    pattern = re.compile(rf"{re.escape(prefix)}_{label}_(\d+)")
    ranked = []
    for name in names:
        match = pattern.fullmatch(str(name))
        rank = int(match.group(1)) if match else 10**9
        ranked.append((rank, str(name)))
    ranked.sort()
    return ranked[0][1]


def choose_canonical_species_name(names: set[str], genus_name: str) -> str:
    # Species numbering is resolved within a given genus namespace.
    if not names:
        raise ValueError("No species names supplied for canonical selection")
    pattern = re.compile(rf"{re.escape(genus_name)} species_(\d+)")
    ranked = []
    for name in names:
        match = pattern.fullmatch(str(name))
        rank = int(match.group(1)) if match else 10**9
        ranked.append((rank, str(name)))
    ranked.sort()
    return ranked[0][1]


def next_custom_index(values: list[str], pattern: re.Pattern[str]) -> int:
    indices = []
    for value in values:
        match = pattern.fullmatch(str(value))
        if match:
            indices.append(int(match.group(1)))
    return max(indices, default=0) + 1


def write_multi_fasta(records: dict[str, object], output_fasta: Path) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([records[genome] for genome in sorted(records)], output_fasta, "fasta")


def run_similarity(combined_records: dict[str, object], output_dir: Path, threads: int) -> Path:
    # A fresh combined taxmyphage similarity run is the authoritative source
    # for deciding whether incoming genomes create new species or genera.
    combined_fasta = output_dir / "combined_genomes.fasta"
    write_multi_fasta(combined_records, combined_fasta)
    command = [
        "taxmyphage",
        "similarity",
        "-i",
        str(combined_fasta),
        "-o",
        str(output_dir),
        "-f",
        "-t",
        str(threads),
        "--no-figures",
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        details = "\n".join(
            part for part in [completed.stdout.strip(), completed.stderr.strip()] if part
        )
        raise RuntimeError(f"taxmyphage similarity failed:\n{details}")
    return output_dir / "similarities.tsv"


def allocate_cluster_names(
    summary_df: pd.DataFrame,
    existing_metadata_df: pd.DataFrame,
    new_genome_ids: set[str],
    prefix: str,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    # Translate raw numeric cluster IDs into stable synthetic names while
    # preserving previous custom names whenever the updated clustering still
    # overlaps them.
    existing_genus_by_genome = existing_metadata_df.set_index("genome")["genus_name"].to_dict()
    existing_species_by_genome = existing_metadata_df.set_index("genome")["species_name"].to_dict()
    existing_family_by_genome = existing_metadata_df.set_index("genome")["family_name"].to_dict()
    existing_representatives = set(
        existing_metadata_df.loc[existing_metadata_df["is_representative"], "genome"]
    )

    next_genus_number = next_custom_index(
        existing_metadata_df["genus_name"].tolist(),
        re.compile(rf"{re.escape(prefix)}_genus_(\d+)"),
    )
    used_species_numbers_by_genus: dict[str, int] = defaultdict(int)
    species_pattern_template = r"{genus} species_(\d+)"
    for _, row in existing_metadata_df.iterrows():
        species_match = re.fullmatch(
            species_pattern_template.format(genus=re.escape(str(row["genus_name"]))),
            str(row["species_name"]),
        )
        if species_match:
            used_species_numbers_by_genus[row["genus_name"]] = max(
                used_species_numbers_by_genus[row["genus_name"]],
                int(species_match.group(1)),
            )

    cluster_df = summary_df.copy()

    genus_name_by_cluster: dict[int, str] = {}
    family_name_by_cluster: dict[int, str] = {}
    species_name_by_cluster: dict[int, str] = {}

    report_rows: list[dict[str, object]] = []

    for genus_cluster, genus_group in cluster_df.groupby("genus_cluster", sort=True):
        genus_genomes = set(genus_group["genome"])
        known_genus_names = {
            existing_genus_by_genome[genome]
            for genome in genus_genomes
            if genome in existing_genus_by_genome
        }
        if known_genus_names:
            genus_name = choose_canonical_name(known_genus_names, prefix, "genus")
            family_names = {
                existing_family_by_genome[genome]
                for genome in genus_genomes
                if genome in existing_family_by_genome
            }
            family_name = next(iter(family_names)) if family_names else family_from_genus_name(genus_name, prefix)
            genus_status = "merged_existing_genera" if len(known_genus_names) > 1 else "existing_genus"
        else:
            genus_name = f"{prefix}_genus_{next_genus_number:03d}"
            family_name = f"{prefix}_family_{next_genus_number:03d}"
            next_genus_number += 1
            genus_status = "new_genus"
        genus_name_by_cluster[genus_cluster] = genus_name
        family_name_by_cluster[genus_cluster] = family_name

        for species_cluster, species_group in genus_group.groupby("species_cluster", sort=True):
            species_genomes = set(species_group["genome"])
            known_species_names = {
                existing_species_by_genome[genome]
                for genome in species_genomes
                if genome in existing_species_by_genome
            }
            # Only species names compatible with the chosen genus name can be
            # preserved. If the genus assignment changes, the species name may
            # need to be reallocated under the new genus namespace.
            canonical_compatible_species_names = {
                species_name
                for species_name in known_species_names
                if str(species_name).startswith(f"{genus_name} ")
            }
            if canonical_compatible_species_names:
                species_name = choose_canonical_species_name(
                    canonical_compatible_species_names, genus_name
                )
                species_status = (
                    "merged_existing_species"
                    if len(canonical_compatible_species_names) > 1 or len(known_species_names) > 1
                    else "existing_species"
                )
            else:
                used_species_numbers_by_genus[genus_name] += 1
                species_name = f"{genus_name} species_{used_species_numbers_by_genus[genus_name]:03d}"
                species_status = "new_species"
            species_name_by_cluster[species_cluster] = species_name

            representative_candidates = species_group.copy()
            existing_rep_rows = representative_candidates[
                representative_candidates["genome"].isin(existing_representatives)
            ]
            # Reuse an existing representative if possible to minimize DB churn.
            if not existing_rep_rows.empty:
                representative_genome = sorted(existing_rep_rows["genome"])[0]
            else:
                representative_genome = (
                    representative_candidates.sort_values(
                        by=["length", "genome"], ascending=[False, True]
                    )
                    .iloc[0]["genome"]
                )

            for row in species_group.itertuples(index=False):
                genome_status = species_status if row.genome in new_genome_ids else "existing_database_member"
                if row.genome in new_genome_ids and genus_status == "new_genus":
                    genome_status = "new_genus"
                report_rows.append(
                    {
                        "genome": row.genome,
                        "genus_cluster": genus_cluster,
                        "species_cluster": species_cluster,
                        "genus_name": genus_name,
                        "species_name": species_name,
                        "family_name": family_name,
                        "status": genome_status,
                        "is_representative": row.genome == representative_genome,
                    }
                )

    assignment_df = pd.DataFrame(report_rows)
    cluster_df = cluster_df.merge(
        assignment_df,
        on=["genome", "genus_cluster", "species_cluster"],
        how="left",
    )
    return cluster_df, report_rows


def build_vmr_from_metadata(
    representative_df: pd.DataFrame,
    existing_vmr_df: pd.DataFrame,
    output_xlsx: Path,
) -> pd.DataFrame:
    # Rebuild the fake VMR workbook to match the updated representative set.
    # Existing rows are reused when available so previous synthetic taxonomy is
    # preserved as much as possible.
    existing_vmr_by_accession = existing_vmr_df.set_index("Virus GENBANK accession").to_dict("index")
    vmr_rows = []
    for isolate_sort, row in enumerate(
        representative_df.sort_values(by=["genus_name", "species_name", "genome"]).itertuples(index=False),
        start=1,
    ):
        if row.genome in existing_vmr_by_accession:
            vmr_row = dict(existing_vmr_by_accession[row.genome])
            vmr_row["Genus"] = row.genus_name
            vmr_row["Species"] = row.species_name
            vmr_row["Family"] = row.family_name
            vmr_row["Virus GENBANK accession"] = row.genome
            vmr_row["Virus name(s)"] = row.genome
            vmr_row["Virus name abbreviation(s)"] = row.genome
            vmr_row["Virus isolate designation"] = row.genome
        else:
            vmr_row = {
                "Isolate ID": f"VMR_CUSTOM_{isolate_sort:06d}",
                "Species Sort": isolate_sort,
                "Isolate Sort": 1,
                "Realm": "CustomRealm",
                "Subrealm": "",
                "Kingdom": "CustomKingdom",
                "Subkingdom": "",
                "Phylum": "CustomPhylum",
                "Subphylum": "",
                "Class": "CustomClass",
                "Subclass": "",
                "Order": "CustomOrder",
                "Suborder": "",
                "Family": row.family_name,
                "Subfamily": "",
                "Genus": row.genus_name,
                "Subgenus": "",
                "Species": row.species_name,
                "ICTV_ID": f"ICTV_CUSTOM_{isolate_sort:06d}",
                "Exemplar or additional isolate": "E",
                "Virus name(s)": row.genome,
                "Virus name abbreviation(s)": row.genome,
                "Virus isolate designation": row.genome,
                "Virus GENBANK accession": row.genome,
                "Genome coverage": "Complete genome",
                "Genome": "dsDNA",
                "Host source": "unknown",
                "Accessions Link": "custom",
            }
        for column in VMR_COLUMNS:
            vmr_row.setdefault(column, "")
        vmr_rows.append({column: vmr_row[column] for column in VMR_COLUMNS})
    vmr_df = pd.DataFrame(vmr_rows, columns=VMR_COLUMNS)
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        vmr_df.to_excel(writer, sheet_name="VMR Custom", index=False)
    return vmr_df


def write_cluster_membership_tables(cluster_df: pd.DataFrame, output_dir: Path) -> None:
    # Persist both user-facing cluster summaries and the canonical all-genome
    # metadata table consumed by future update runs.
    cluster_df = cluster_df.sort_values(
        by=["genus_name", "species_name", "genome"]
    ).reset_index(drop=True)
    cluster_df.to_csv(output_dir / "cluster_assignments.tsv", sep="\t", index=False)
    cluster_df.to_csv(output_dir / "all_genomes_metadata.tsv", sep="\t", index=False)

    genus_df = (
        cluster_df.groupby(["genus_cluster", "genus_name", "family_name"], sort=True)["genome"]
        .apply(lambda values: ",".join(values))
        .reset_index(name="genomes")
    )
    genus_df.to_csv(output_dir / "genus_clusters.tsv", sep="\t", index=False)

    species_df = (
        cluster_df.groupby(
            ["genus_cluster", "genus_name", "species_cluster", "species_name", "family_name"],
            sort=True,
        )["genome"]
        .apply(lambda values: ",".join(values))
        .reset_index(name="genomes")
    )
    species_df.to_csv(output_dir / "species_clusters.tsv", sep="\t", index=False)


def classify_new_genome_statuses(cluster_df: pd.DataFrame, new_genome_ids: set[str]) -> pd.DataFrame:
    # This report is the main per-input summary returned by the updater.
    new_df = cluster_df[cluster_df["genome"].isin(new_genome_ids)].copy()
    if new_df.empty:
        return new_df
    new_df = new_df[
        [
            "genome",
            "genus_cluster",
            "species_cluster",
            "genus_name",
            "species_name",
            "family_name",
            "status",
            "is_representative",
            "length",
        ]
    ].sort_values(by=["status", "genus_name", "species_name", "genome"])
    return new_df


def main() -> int:
    args = parse_args()

    db_dir = args.db_dir.resolve()
    output_db_dir = (
        args.output_db_dir.resolve() if args.output_db_dir else db_dir
    )
    input_fastas = find_input_fastas(args.input)

    existing_records = load_existing_db_records(db_dir)
    new_records = load_records_from_fastas(input_fastas)

    overlapping_ids = set(existing_records).intersection(new_records)
    if overlapping_ids:
        for genome_id in sorted(overlapping_ids):
            print(f"Skipping existing genome already present in the database: {genome_id}")
            new_records.pop(genome_id, None)

    if not new_records:
        print("No new genomes remain after skipping inputs already present in the database.")
        return 0

    existing_metadata_df, existing_vmr_df = load_existing_metadata(db_dir, args.prefix)
    combined_records = dict(existing_records)
    combined_records.update(new_records)

    # Similarity is computed in a temporary workspace first. The real database
    # directory is only updated after the new clustering has been fully derived.
    temp_dir_handle = tempfile.TemporaryDirectory(prefix="custom_db_update_")
    work_dir = Path(temp_dir_handle.name)
    try:
        run_similarity(combined_records, work_dir, args.threads)
        similarities_df = pd.read_csv(work_dir / "similarities.tsv", sep="\t")
        summary_df, _, _ = build_summary_tables(
            records=combined_records,
            similarities_df=similarities_df,
            genus_threshold=70.0,
            species_threshold=95.0,
        )
        cluster_df, _ = allocate_cluster_names(
            summary_df=summary_df,
            existing_metadata_df=existing_metadata_df,
            new_genome_ids=set(new_records),
            prefix=args.prefix,
        )

        new_genome_report_df = classify_new_genome_statuses(cluster_df, set(new_records))
        updated_representatives_df = cluster_df[cluster_df["is_representative"]].copy()

        if output_db_dir != db_dir and output_db_dir.exists():
            raise FileExistsError(
                f"Output database directory already exists: {output_db_dir}. "
                "Use a fresh directory or omit --output-db-dir to update in place."
            )
        output_db_dir.mkdir(parents=True, exist_ok=True)

        all_genomes_fasta = output_db_dir / "all_genomes.fasta"
        representatives_fasta = output_db_dir / "representatives.fasta"
        db_fasta = output_db_dir / "Bacteriophage_genomes.fasta"
        db_fasta_gz = output_db_dir / "Bacteriophage_genomes.fasta.gz"
        vmr_path = output_db_dir / "VMR.xlsx"
        msh_path = output_db_dir / "ICTV.msh"
        mpa_path = output_db_dir / "M.pa"

        write_multi_fasta(combined_records, all_genomes_fasta)
        # Keep all genomes for future update cycles, but still build the actual
        # queryable taxmyphage reference DB from one representative per species.
        write_representatives_fasta(updated_representatives_df, combined_records, representatives_fasta)
        shutil.copyfile(representatives_fasta, db_fasta)
        gzip_copy(db_fasta, db_fasta_gz)
        make_blast_db(db_fasta)

        vmr_df = build_vmr_from_metadata(
            representative_df=updated_representatives_df,
            existing_vmr_df=existing_vmr_df,
            output_xlsx=vmr_path,
        )
        mash_reference_fastas = stage_mash_references(
            representatives_df=updated_representatives_df,
            vmr_df=vmr_df,
            records=combined_records,
            mash_reference_dir=output_db_dir / "mash_references",
        )
        make_mash_sketch(mash_reference_fastas, msh_path, args.threads)
        if args.build_mpa:
            idab_df = build_precomputed_blast_parquet(representatives_fasta, mpa_path, args.threads)
            write_pairwise_similarity_table(
                idab_df=idab_df,
                records={row.genome: combined_records[row.genome] for row in updated_representatives_df.itertuples(index=False)},
                output_tsv=output_db_dir / "Mpa_pairwise_similarity.tsv",
            )
        else:
            if mpa_path.exists():
                mpa_path.unlink()
            pairwise_tsv = output_db_dir / "Mpa_pairwise_similarity.tsv"
            if pairwise_tsv.exists():
                pairwise_tsv.unlink()

        write_cluster_membership_tables(cluster_df, output_db_dir)
        updated_representatives_df.to_csv(output_db_dir / "representatives.tsv", sep="\t", index=False)
        new_genome_report_df.to_csv(output_db_dir / "new_genome_report.tsv", sep="\t", index=False)
        similarities_df.to_csv(output_db_dir / "last_update_similarities.tsv", sep="\t", index=False)
        write_readme(output_db_dir / "README.txt")

        num_new_species = int((new_genome_report_df["status"] == "new_species").sum())
        num_new_genera = int((new_genome_report_df["status"] == "new_genus").sum())
        print(f"Updated custom database written to {output_db_dir}")
        print(f"Input genomes: {len(new_records)}")
        print(f"New species detected: {num_new_species}")
        print(f"New genera detected: {num_new_genera}")
        print("Use --no-precomputed when querying this custom database.")
        if not args.build_mpa:
            print("M.pa was not built, so taxmyphage will default away from precomputed mode.")
    finally:
        if args.keep_workdir:
            print(f"Kept temporary workdir at {work_dir}")
            temp_dir_handle.cleanup = lambda: None
        else:
            temp_dir_handle.cleanup()

    return 0


if __name__ == "__main__":
    sys.exit(main())
