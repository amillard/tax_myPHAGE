#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# This script builds a self-contained custom taxmyphage database from:
# 1. local FASTA files
# 2. an existing taxmyphage similarity run containing similarities.tsv
#
# The core idea is:
# - cluster genomes at genus/species thresholds
# - choose one representative per species cluster
# - rebuild the assets taxmyphage expects (FASTA, BLAST DB, MASH sketch, VMR, M.pa)
#
# M.pa stores ordered-pair identity counts idAB among database representatives.
# The final pairwise similarity is derived later using:
#   simAB = ((idAB + idBA) * 100) / (lA + lB)
#   distAB = 100 - simAB
#
# A validation table with these explicit calculations is also written out.

GENUS_THRESHOLD = 70.0
SPECIES_THRESHOLD = 95.0
FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fsa"}
VMR_COLUMNS = [
    "Isolate ID",
    "Species Sort",
    "Isolate Sort",
    "Realm",
    "Subrealm",
    "Kingdom",
    "Subkingdom",
    "Phylum",
    "Subphylum",
    "Class",
    "Subclass",
    "Order",
    "Suborder",
    "Family",
    "Subfamily",
    "Genus",
    "Subgenus",
    "Species",
    "ICTV_ID",
    "Exemplar or additional isolate",
    "Virus name(s)",
    "Virus name abbreviation(s)",
    "Virus isolate designation",
    "Virus GENBANK accession",
    "Genome coverage",
    "Genome",
    "Host source",
    "Accessions Link",
]


class UnionFind:
    # Minimal connected-component helper used for threshold clustering.
    # A genome graph is built implicitly by unioning every pair with sim >= threshold.
    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}
        self.rank = {item: 0 for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1

    def clusters(self) -> list[list[str]]:
        members: dict[str, list[str]] = defaultdict(list)
        for item in self.parent:
            members[self.find(item)].append(item)
        return sorted((sorted(cluster) for cluster in members.values()), key=lambda c: (c[0], len(c)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a custom taxmyphage database from local FASTA files and "
            "a taxmyphage similarity output folder."
        )
    )
    parser.add_argument(
        "--fasta-dir",
        type=Path,
        default=Path("."),
        help="Directory containing phage FASTA files. Default: current directory.",
    )
    parser.add_argument(
        "--similarity-dir",
        type=Path,
        default=Path("taxmyphage_similarity"),
        help="Directory containing similarities.tsv. Default: taxmyphage_similarity.",
    )
    parser.add_argument(
        "--db-dir",
        type=Path,
        default=Path("custom_taxmyphage_db"),
        help="Output directory for the custom taxmyphage database.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for BLAST and mash sketch. Default: 1.",
    )
    parser.add_argument(
        "--species-threshold",
        type=float,
        default=SPECIES_THRESHOLD,
        help="Species threshold in percent similarity. Default: 95.",
    )
    parser.add_argument(
        "--genus-threshold",
        type=float,
        default=GENUS_THRESHOLD,
        help="Genus threshold in percent similarity. Default: 70.",
    )
    parser.add_argument(
        "--prefix",
        default="Custom",
        help="Prefix for fake taxon names in the custom VMR. Default: Custom.",
    )
    parser.add_argument(
        "--build-mpa",
        action="store_true",
        help="Build M.pa and Mpa_pairwise_similarity.tsv. Default: off for safer no-precomputed querying.",
    )
    return parser.parse_args()


def find_fasta_files(fasta_dir: Path) -> list[Path]:
    # The build script assumes one genome per file because the genome ID is
    # taken from the filename stem and reused throughout the custom DB.
    fasta_files = [
        path
        for path in sorted(fasta_dir.iterdir())
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    ]
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {fasta_dir}")
    return fasta_files


def load_records(fasta_files: list[Path]) -> dict[str, object]:
    # Normalize each FASTA record so the sequence ID used everywhere in the DB
    # is the filename stem rather than any original FASTA header variant.
    records = {}
    for fasta_path in fasta_files:
        entries = list(SeqIO.parse(str(fasta_path), "fasta"))
        if len(entries) != 1:
            raise ValueError(
                f"{fasta_path} contains {len(entries)} FASTA records. "
                "This script expects one genome per FASTA file."
            )
        record = entries[0]
        genome_id = fasta_path.stem
        record.id = genome_id
        record.name = genome_id
        record.description = genome_id
        records[genome_id] = record
    return records


def load_similarity_table(similarity_dir: Path) -> pd.DataFrame:
    similarities_path = similarity_dir / "similarities.tsv"
    if not similarities_path.exists():
        raise FileNotFoundError(f"Could not find {similarities_path}")
    df = pd.read_csv(similarities_path, sep="\t")
    missing = {"A", "B", "sim"} - set(df.columns)
    if missing:
        raise ValueError(f"{similarities_path} is missing columns: {sorted(missing)}")
    return df


def build_cluster_assignments(
    genomes: list[str],
    similarities_df: pd.DataFrame,
    threshold: float,
) -> dict[str, int]:
    # Genus/species are both produced as connected components in the pairwise
    # similarity graph. This means transitive relationships are respected:
    # if A~B and B~C at a threshold, they become one cluster.
    uf = UnionFind(genomes)
    for row in similarities_df.itertuples(index=False):
        if row.sim >= threshold:
            uf.union(row.A, row.B)
    assignments = {}
    for cluster_id, members in enumerate(uf.clusters(), start=1):
        for genome in members:
            assignments[genome] = cluster_id
    return assignments


def choose_representatives(
    summary_df: pd.DataFrame, records: dict[str, object]
) -> pd.DataFrame:
    # One representative is chosen per species cluster.
    # Current rule: longest genome wins, then genome name as a stable tie-breaker.
    picked_rows = []
    for _, species_df in summary_df.groupby("species_cluster", sort=True):
        species_df = species_df.copy()
        species_df["length"] = species_df["genome"].map(lambda genome: len(records[genome].seq))
        species_df = species_df.sort_values(
            by=["length", "genome"], ascending=[False, True]
        ).reset_index(drop=True)
        chosen = species_df.iloc[0].to_dict()
        chosen["representative"] = True
        picked_rows.append(chosen)
    return pd.DataFrame(picked_rows).sort_values(
        by=["genus_cluster", "species_cluster", "genome"]
    ).reset_index(drop=True)


def build_summary_tables(
    records: dict[str, object],
    similarities_df: pd.DataFrame,
    genus_threshold: float,
    species_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Build the full per-genome assignment table, then derived genus/species
    # membership tables that are easier to inspect downstream.
    genomes = sorted(records)
    genus_assignment = build_cluster_assignments(genomes, similarities_df, genus_threshold)
    species_assignment = build_cluster_assignments(genomes, similarities_df, species_threshold)

    rows = []
    for genome in genomes:
        rows.append(
            {
                "genome": genome,
                "length": len(records[genome].seq),
                "genus_cluster": genus_assignment[genome],
                "species_cluster": species_assignment[genome],
            }
        )
    summary_df = pd.DataFrame(rows).sort_values(
        by=["genus_cluster", "species_cluster", "genome"]
    ).reset_index(drop=True)

    genus_members = (
        summary_df.groupby("genus_cluster")["genome"]
        .apply(lambda genomes_in_cluster: ",".join(genomes_in_cluster))
        .reset_index(name="genomes")
    )
    species_members = (
        summary_df.groupby(["genus_cluster", "species_cluster"])["genome"]
        .apply(lambda genomes_in_cluster: ",".join(genomes_in_cluster))
        .reset_index(name="genomes")
    )
    return summary_df, genus_members, species_members


def write_representatives_fasta(
    representatives_df: pd.DataFrame,
    records: dict[str, object],
    output_fasta: Path,
) -> None:
    # The representative FASTA is the actual sequence set used to construct
    # the custom BLAST DB, MASH sketch, M.pa, and VMR.
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    selected_records = [records[row.genome] for row in representatives_df.itertuples(index=False)]
    SeqIO.write(selected_records, output_fasta, "fasta")


def gzip_copy(source: Path, destination: Path) -> None:
    with open(source, "rb") as src, gzip.open(destination, "wb") as dst:
        shutil.copyfileobj(src, dst)


def run_command(command: list[str], description: str) -> None:
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        stdout = completed.stdout.strip()
        details = "\n".join(part for part in [stdout, stderr] if part)
        raise RuntimeError(f"{description} failed:\n{details}")


def make_blast_db(db_fasta: Path) -> None:
    run_command(
        [
            "makeblastdb",
            "-in",
            str(db_fasta),
            "-parse_seqids",
            "-dbtype",
            "nucl",
        ],
        "makeblastdb",
    )


def build_precomputed_blast_parquet(
    representative_fasta: Path,
    output_parquet: Path,
    threads: int,
) -> pd.DataFrame:
    # Build M.pa from an all-vs-all BLAST of representative genomes.
    #
    # Important implementation detail:
    # idAB is computed in query-space for ordered pair A -> B.
    # We collect identical aligned query positions and store each query position
    # only once, which collapses overlap among multiple BLAST HSPs.
    #
    # This is the cached quantity taxmyphage later combines with idBA and
    # genome lengths to calculate pairwise similarity.
    with tempfile.TemporaryDirectory(prefix="custom_taxmyphage_") as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        blast_db = temp_dir / "repdb"
        blast_tsv = temp_dir / "representatives.blast.tsv"

        run_command(
            [
                "makeblastdb",
                "-in",
                str(representative_fasta),
                "-out",
                str(blast_db),
                "-dbtype",
                "nucl",
            ],
            "makeblastdb for M.pa",
        )

        blast_cmd = [
            "blastn",
            "-evalue",
            "1",
            "-max_target_seqs",
            "10000",
            "-num_threads",
            str(threads),
            "-word_size",
            "7",
            "-reward",
            "2",
            "-penalty",
            "-3",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-query",
            str(representative_fasta),
            "-db",
            str(blast_db),
            "-outfmt",
            "6 qseqid sseqid pident length qlen slen mismatch nident gapopen qstart qend sstart send qseq sseq evalue bitscore",
            "-out",
            str(blast_tsv),
        ]
        run_command(blast_cmd, "blastn for M.pa")

        identity_positions_by_pair: dict[tuple[str, str], set[int]] = defaultdict(set)

        with open(blast_tsv, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) != 17:
                    raise ValueError(
                        f"Unexpected BLAST output format in {blast_tsv}: {raw_line[:120]}"
                    )
                (
                    qseqid,
                    sseqid,
                    _pident,
                    _length,
                    _qlen,
                    _slen,
                    _mismatch,
                    _nident,
                    _gapopen,
                    qstart,
                    qend,
                    _sstart,
                    _send,
                    qseq,
                    sseq,
                    _evalue,
                    _bitscore,
                ) = fields
                pair = (qseqid, sseqid)
                start = int(qstart)
                end = int(qend)
                step = 1 if end >= start else -1
                # Positions are tracked in query coordinates. That makes this
                # ordered-pair specific and is the reason idAB and idBA are
                # similar but not guaranteed to be identical.
                for position, (query_base, subject_base) in zip(
                    range(start, end + step, step), zip(qseq, sseq)
                ):
                    if query_base != "-" and query_base == subject_base:
                        identity_positions_by_pair[pair].add(position)

        if not identity_positions_by_pair:
            raise RuntimeError("No BLAST identities were parsed while building M.pa")

        parquet_df = pd.DataFrame(
            (
                {"A": genome_a, "B": genome_b, "idAB": len(identity_positions)}
                for (genome_a, genome_b), identity_positions in sorted(identity_positions_by_pair.items())
            )
        )
        parquet_df.to_parquet(output_parquet, index=False)
        return parquet_df


def write_pairwise_similarity_table(
    idab_df: pd.DataFrame,
    records: dict[str, object],
    output_tsv: Path,
) -> None:
    # This file is written for transparency and validation.
    # taxmyphage consumes M.pa, but this table shows the explicit formula:
    # simAB = ((idAB + idBA) * 100) / (lA + lB)
    # distAB = 100 - simAB
    idab_lookup = idab_df.set_index(["A", "B"])["idAB"].to_dict()
    genomes = sorted(records)
    rows = []
    for i, genome_a in enumerate(genomes):
        for genome_b in genomes[i:]:
            id_ab = int(idab_lookup.get((genome_a, genome_b), 0))
            id_ba = int(idab_lookup.get((genome_b, genome_a), 0))
            len_a = len(records[genome_a].seq)
            len_b = len(records[genome_b].seq)
            sim_ab = ((id_ab + id_ba) * 100.0) / (len_a + len_b)
            dist_ab = 100.0 - sim_ab
            rows.append(
                {
                    "A": genome_a,
                    "B": genome_b,
                    "idAB": id_ab,
                    "idBA": id_ba,
                    "lA": len_a,
                    "lB": len_b,
                    "simAB": sim_ab,
                    "distAB": dist_ab,
                }
            )
    pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)


def build_fake_vmr(
    representatives_df: pd.DataFrame,
    output_xlsx: Path,
    prefix: str,
) -> pd.DataFrame:
    # Build the minimum workbook taxmyphage needs for a custom reference set.
    # Names are synthetic placeholders; the key matching field is the accession
    # column, which must match the FASTA / BLAST / MASH identifiers.
    genus_name_map = {
        cluster_id: f"{prefix}_genus_{cluster_id:03d}"
        for cluster_id in sorted(representatives_df["genus_cluster"].unique())
    }

    species_order_within_genus: dict[int, dict[int, int]] = defaultdict(dict)
    for genus_cluster, species_cluster in (
        representatives_df[["genus_cluster", "species_cluster"]]
        .drop_duplicates()
        .sort_values(by=["genus_cluster", "species_cluster"])
        .itertuples(index=False)
    ):
        species_order_within_genus[genus_cluster][species_cluster] = (
            len(species_order_within_genus[genus_cluster]) + 1
        )

    vmr_rows = []
    for isolate_sort, row in enumerate(
        representatives_df.sort_values(by=["genus_cluster", "species_cluster", "genome"]).itertuples(index=False),
        start=1,
    ):
        genus_name = genus_name_map[row.genus_cluster]
        species_index = species_order_within_genus[row.genus_cluster][row.species_cluster]
        species_name = f"{genus_name} species_{species_index:03d}"
        vmr_rows.append(
            {
                "Isolate ID": f"VMR_CUSTOM_{isolate_sort:06d}",
                "Species Sort": species_index,
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
                "Family": f"{prefix}_family_{row.genus_cluster:03d}",
                "Subfamily": "",
                "Genus": genus_name,
                "Subgenus": "",
                "Species": species_name,
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
        )

    vmr_df = pd.DataFrame(vmr_rows, columns=VMR_COLUMNS)
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        vmr_df.to_excel(writer, sheet_name="VMR Custom", index=False)
    return vmr_df


def write_readme(output_path: Path) -> None:
    output_path.write_text(
        "\n".join(
            [
                "Custom taxmyphage database",
                "",
                "Important:",
                "Run taxmyphage against this custom database with --no-precomputed.",
                "",
                "Reason:",
                "taxmyphage precomputed mode can overestimate similarity for a new query",
                "when reverse query/reference identities are unavailable in the cached matrix.",
                "",
                "Example:",
                "taxmyphage run -i QUERY.fa -o results -t 12 --no-figures --no-precomputed -db custom_taxmyphage_db",
                "",
            ]
        ),
        encoding="utf-8",
    )


def stage_mash_references(
    representatives_df: pd.DataFrame,
    vmr_df: pd.DataFrame,
    records: dict[str, object],
    mash_reference_dir: Path,
) -> list[Path]:
    # taxmyphage expects MASH references to look like per-genome files in
    # genus-specific directories. Using a single combined FASTA would lose the
    # accession-level identity needed to map MASH hits back to the VMR rows.
    mash_reference_dir.mkdir(parents=True, exist_ok=True)
    genus_by_accession = vmr_df.set_index("Virus GENBANK accession")["Genus"].to_dict()
    staged_files = []
    for row in representatives_df.sort_values(by=["genus_cluster", "species_cluster", "genome"]).itertuples(index=False):
        genus_name = genus_by_accession[row.genome]
        genus_dir = mash_reference_dir / genus_name
        genus_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = genus_dir / f"{row.genome}.fasta"
        SeqIO.write([records[row.genome]], fasta_path, "fasta")
        staged_files.append(fasta_path)
    return staged_files


def make_mash_sketch(reference_fastas: list[Path], msh_path: Path, threads: int) -> None:
    if not reference_fastas:
        raise ValueError("No representative FASTA files were staged for mash sketch")
    run_command(
        [
            "mash",
            "sketch",
            "-p",
            str(threads),
            "-o",
            str(msh_path.with_suffix("")),
            *[str(path) for path in reference_fastas],
        ],
        "mash sketch",
    )


def main() -> int:
    args = parse_args()

    fasta_dir = args.fasta_dir.resolve()
    similarity_dir = args.similarity_dir.resolve()
    db_dir = args.db_dir.resolve()
    db_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = find_fasta_files(fasta_dir)
    records = load_records(fasta_files)
    similarities_df = load_similarity_table(similarity_dir)

    summary_df, genus_members_df, species_members_df = build_summary_tables(
        records=records,
        similarities_df=similarities_df,
        genus_threshold=args.genus_threshold,
        species_threshold=args.species_threshold,
    )
    representatives_df = choose_representatives(summary_df, records)

    # These tables document the clustering state used to define the custom DB.
    summary_df.to_csv(db_dir / "cluster_assignments.tsv", sep="\t", index=False)
    genus_members_df.to_csv(db_dir / "genus_clusters.tsv", sep="\t", index=False)
    species_members_df.to_csv(db_dir / "species_clusters.tsv", sep="\t", index=False)
    representatives_df.to_csv(db_dir / "representatives.tsv", sep="\t", index=False)

    db_fasta = db_dir / "Bacteriophage_genomes.fasta"
    db_fasta_gz = db_dir / "Bacteriophage_genomes.fasta.gz"
    representative_fasta = db_dir / "representatives.fasta"
    msh_path = db_dir / "ICTV.msh"
    vmr_path = db_dir / "VMR.xlsx"
    precomputed_path = db_dir / "M.pa"

    write_representatives_fasta(representatives_df, records, representative_fasta)
    # taxmyphage expects the BLAST DB to be called Bacteriophage_genomes.fasta.
    shutil.copyfile(representative_fasta, db_fasta)
    gzip_copy(db_fasta, db_fasta_gz)
    make_blast_db(db_fasta)
    vmr_df = build_fake_vmr(representatives_df, vmr_path, args.prefix)
    mash_reference_fastas = stage_mash_references(
        representatives_df=representatives_df,
        vmr_df=vmr_df,
        records=records,
        mash_reference_dir=db_dir / "mash_references",
    )
    make_mash_sketch(mash_reference_fastas, msh_path, args.threads)
    if args.build_mpa:
        idab_df = build_precomputed_blast_parquet(representative_fasta, precomputed_path, args.threads)
        write_pairwise_similarity_table(
            idab_df=idab_df,
            records={row.genome: records[row.genome] for row in representatives_df.itertuples(index=False)},
            output_tsv=db_dir / "Mpa_pairwise_similarity.tsv",
        )
    else:
        if precomputed_path.exists():
            precomputed_path.unlink()
        pairwise_tsv = db_dir / "Mpa_pairwise_similarity.tsv"
        if pairwise_tsv.exists():
            pairwise_tsv.unlink()
    write_readme(db_dir / "README.txt")

    print(f"Wrote custom taxmyphage database to {db_dir}")
    print(f"Representatives: {len(representatives_df)}")
    print(f"Genera: {summary_df['genus_cluster'].nunique()}")
    print(f"Species: {summary_df['species_cluster'].nunique()}")
    print("Use --no-precomputed when querying this custom database.")
    if not args.build_mpa:
        print("M.pa was not built, so taxmyphage will default away from precomputed mode.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
