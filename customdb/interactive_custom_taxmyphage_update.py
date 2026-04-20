#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import pandas as pd

from build_custom_taxmyphage_db import (
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
from update_custom_taxmyphage_db import (
    allocate_cluster_names,
    build_vmr_from_metadata,
    classify_new_genome_statuses,
    find_input_fastas,
    load_existing_db_records,
    load_existing_metadata,
    load_records_from_fastas,
    run_similarity,
    write_cluster_membership_tables,
    write_multi_fasta,
)


# Interactive wrapper for maintaining a custom taxmyphage database.
#
# The script first builds a complete candidate updated database outside the
# live custom DB directory. The live DB is only replaced after the user has
# reviewed the new-genome report and confirmed the update.

GENERATED_ROOT_FILES = {
    "Bacteriophage_genomes.fasta",
    "Bacteriophage_genomes.fasta.gz",
    "ICTV.msh",
    "M.pa",
    "Mpa_pairwise_similarity.tsv",
    "README.txt",
    "VMR.xlsx",
    "all_genomes.fasta",
    "all_genomes_metadata.tsv",
    "cluster_assignments.tsv",
    "genus_clusters.tsv",
    "last_update_similarities.tsv",
    "new_genome_report.tsv",
    "representatives.fasta",
    "representatives.tsv",
    "species_clusters.tsv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Interactively compare new genomes against a custom taxmyphage DB, "
            "then apply VMR and database updates after confirmation."
        )
    )
    parser.add_argument(
        "--db-dir",
        type=Path,
        default=Path("custom_taxmyphage_db"),
        help="Existing custom database directory. Default: custom_taxmyphage_db.",
    )
    parser.add_argument(
        "--input",
        nargs="*",
        default=None,
        help="Input FASTA files or directories. If omitted, prompt interactively.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for taxmyphage similarity, BLAST, and mash sketch. Default: 1.",
    )
    parser.add_argument(
        "--prefix",
        default="Custom",
        help="Prefix used for synthetic family/genus/species names. Default: Custom.",
    )
    parser.add_argument(
        "--candidate-db-dir",
        type=Path,
        default=None,
        help="Optional directory for the candidate updated DB. Must not exist.",
    )
    parser.add_argument(
        "--keep-candidate",
        action="store_true",
        help="Keep the candidate DB directory instead of deleting it after applying/skipping.",
    )
    parser.add_argument(
        "--keep-workdir",
        action="store_true",
        help="Keep the temporary taxmyphage similarity work directory for inspection.",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Apply without interactive confirmation when new species or genera are found.",
    )
    parser.add_argument(
        "--apply-existing-species",
        action="store_true",
        help=(
            "Also apply updates when inputs are new genomes in existing species. "
            "By default, only new species or genera are auto-eligible."
        ),
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not create a timestamped backup of the existing DB before applying.",
    )
    parser.add_argument(
        "--build-mpa",
        action="store_true",
        help="Build M.pa and Mpa_pairwise_similarity.tsv. Default: off.",
    )
    return parser.parse_args()


def prompt_for_inputs() -> list[str]:
    print("Enter FASTA files or directories to check.")
    print("Submit a blank line when done.")
    inputs: list[str] = []
    while True:
        raw_value = input("Input path: ").strip()
        if not raw_value:
            break
        inputs.append(raw_value)
    return inputs


def prompt_yes_no(question: str, default: bool = False) -> bool:
    suffix = " [Y/n]: " if default else " [y/N]: "
    while True:
        answer = input(question + suffix).strip().lower()
        if not answer:
            return default
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please answer yes or no.")


def create_candidate_dir(args: argparse.Namespace) -> Path:
    if args.candidate_db_dir is not None:
        candidate_dir = args.candidate_db_dir.resolve()
        if candidate_dir.exists():
            raise FileExistsError(f"Candidate DB directory already exists: {candidate_dir}")
        candidate_dir.mkdir(parents=True)
        return candidate_dir

    return Path(tempfile.mkdtemp(prefix="custom_taxmyphage_candidate_"))


def build_candidate_database(
    db_dir: Path,
    input_fastas: list[Path],
    candidate_dir: Path,
    threads: int,
    prefix: str,
    build_mpa: bool,
    keep_workdir: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, Path | None]:
    existing_records = load_existing_db_records(db_dir)
    new_records = load_records_from_fastas(input_fastas)

    overlapping_ids = set(existing_records).intersection(new_records)
    for genome_id in sorted(overlapping_ids):
        print(f"Skipping existing genome already present in the database: {genome_id}")
        new_records.pop(genome_id, None)

    if not new_records:
        return pd.DataFrame(), pd.DataFrame(), None

    existing_metadata_df, existing_vmr_df = load_existing_metadata(db_dir, prefix)
    combined_records = dict(existing_records)
    combined_records.update(new_records)

    work_dir = Path(tempfile.mkdtemp(prefix="custom_db_update_similarity_"))
    kept_work_dir: Path | None = None
    try:
        run_similarity(combined_records, work_dir, threads)
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
            prefix=prefix,
        )

        new_genome_report_df = classify_new_genome_statuses(cluster_df, set(new_records))
        updated_representatives_df = cluster_df[cluster_df["is_representative"]].copy()

        all_genomes_fasta = candidate_dir / "all_genomes.fasta"
        representatives_fasta = candidate_dir / "representatives.fasta"
        db_fasta = candidate_dir / "Bacteriophage_genomes.fasta"
        db_fasta_gz = candidate_dir / "Bacteriophage_genomes.fasta.gz"
        vmr_path = candidate_dir / "VMR.xlsx"
        msh_path = candidate_dir / "ICTV.msh"
        mpa_path = candidate_dir / "M.pa"

        write_multi_fasta(combined_records, all_genomes_fasta)
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
            mash_reference_dir=candidate_dir / "mash_references",
        )
        make_mash_sketch(mash_reference_fastas, msh_path, threads)

        if build_mpa:
            idab_df = build_precomputed_blast_parquet(representatives_fasta, mpa_path, threads)
            representative_records = {
                row.genome: combined_records[row.genome]
                for row in updated_representatives_df.itertuples(index=False)
            }
            write_pairwise_similarity_table(
                idab_df=idab_df,
                records=representative_records,
                output_tsv=candidate_dir / "Mpa_pairwise_similarity.tsv",
            )

        write_cluster_membership_tables(cluster_df, candidate_dir)
        updated_representatives_df.to_csv(candidate_dir / "representatives.tsv", sep="\t", index=False)
        new_genome_report_df.to_csv(candidate_dir / "new_genome_report.tsv", sep="\t", index=False)
        similarities_df.to_csv(candidate_dir / "last_update_similarities.tsv", sep="\t", index=False)
        write_readme(candidate_dir / "README.txt")

        if keep_workdir:
            kept_work_dir = work_dir
    finally:
        if not keep_workdir and work_dir.exists():
            shutil.rmtree(work_dir)

    return new_genome_report_df, cluster_df, kept_work_dir


def display_report(report_df: pd.DataFrame) -> None:
    if report_df.empty:
        print("No new genomes remain after skipping inputs already present in the database.")
        return

    columns = [
        "genome",
        "status",
        "genus_name",
        "species_name",
        "is_representative",
        "length",
    ]
    print("\nNew-genome report:")
    print(report_df[columns].to_string(index=False))


def should_apply_update(report_df: pd.DataFrame, args: argparse.Namespace) -> bool:
    statuses = set(report_df["status"]) if not report_df.empty else set()
    has_new_taxon = bool(statuses.intersection({"new_genus", "new_species"}))

    if args.yes:
        return has_new_taxon or args.apply_existing_species

    if has_new_taxon:
        return prompt_yes_no("Apply this update to the existing custom database?", default=True)

    if args.apply_existing_species:
        return prompt_yes_no(
            "No new species/genera were found, but new genomes can be recorded. Apply anyway?",
            default=False,
        )

    print("\nNo new genus or species was detected, so the live custom DB was not changed.")
    print("Use --apply-existing-species if you also want to record same-species genomes.")
    return False


def backup_database(db_dir: Path) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = db_dir.with_name(f"{db_dir.name}.backup_{timestamp}")
    shutil.copytree(db_dir, backup_dir)
    return backup_dir


def clear_generated_outputs(db_dir: Path) -> None:
    for filename in GENERATED_ROOT_FILES:
        path = db_dir / filename
        if path.exists():
            path.unlink()
    for path in db_dir.glob("Bacteriophage_genomes.fasta.n*"):
        if path.is_file():
            path.unlink()
    mash_references = db_dir / "mash_references"
    if mash_references.exists():
        shutil.rmtree(mash_references)


def copy_candidate_to_database(candidate_dir: Path, db_dir: Path) -> None:
    clear_generated_outputs(db_dir)
    for source in candidate_dir.iterdir():
        destination = db_dir / source.name
        if source.is_dir():
            shutil.copytree(source, destination, dirs_exist_ok=True)
        else:
            shutil.copy2(source, destination)


def main() -> int:
    args = parse_args()
    db_dir = args.db_dir.resolve()
    if not db_dir.exists():
        raise FileNotFoundError(f"Custom database directory does not exist: {db_dir}")

    raw_inputs = args.input if args.input else prompt_for_inputs()
    if not raw_inputs:
        print("No input FASTA files or directories were supplied.")
        return 1

    input_fastas = find_input_fastas(raw_inputs)
    print(f"Checking {len(input_fastas)} input FASTA file(s) against {db_dir}")

    candidate_dir = create_candidate_dir(args)
    remove_candidate = not args.keep_candidate
    try:
        report_df, _cluster_df, kept_work_dir = build_candidate_database(
            db_dir=db_dir,
            input_fastas=input_fastas,
            candidate_dir=candidate_dir,
            threads=args.threads,
            prefix=args.prefix,
            build_mpa=args.build_mpa,
            keep_workdir=args.keep_workdir,
        )
        display_report(report_df)

        if kept_work_dir is not None:
            print(f"Kept similarity work directory: {kept_work_dir}")

        if report_df.empty:
            return 0

        if not should_apply_update(report_df, args):
            print(f"Candidate database retained at: {candidate_dir}")
            remove_candidate = False
            return 0

        backup_dir = None
        if not args.no_backup:
            backup_dir = backup_database(db_dir)
            print(f"Backup written to: {backup_dir}")

        copy_candidate_to_database(candidate_dir, db_dir)
        print(f"Updated custom database written to: {db_dir}")
        print("Use --no-precomputed when querying this custom database.")
        return 0
    finally:
        if remove_candidate and candidate_dir.exists():
            shutil.rmtree(candidate_dir)


if __name__ == "__main__":
    sys.exit(main())
