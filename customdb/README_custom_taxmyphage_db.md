# Custom Database Workflow


Two scritpts to build a custom `taxmyphage` database built from local phage genomes.
With a wrapper script 
Still experimental!! warning !! 

## Files

- `build_custom_taxmyphage_db.py`
  Builds a custom database from a set of FASTA files plus a `taxmyphage similarity` result folder containing `similarities.tsv`

- Requires `taxmyphage similiartiy` to be run on the set of genomes first 

- `update_custom_taxmyphage_db.py`
  Takes an existing custom database and then compares against the existing custom database file. Identifies a new generea/species and rebuilds the database updating the mash, database and fake VMR file 
  
  

## Core rules used

- Same genus: similarity `>= 70%`
- Same species: similarity `>= 95%`

These thresholds are applied as connected-component clustering on the pairwise similarity graph:

- edges with `sim >= 70` define genus clusters
- edges with `sim >= 95` define species clusters

Can add in any values to get strains if required. Just set the values 



Use custom databases created by these scripts with:

```bash
taxmyphage run -i QUERY.fa -o results --no-figures --no-precomputed -db custom_taxmyphage_db
```

`M.pa` is not created correctly currently 



##### What `build_custom_taxmyphage_db.py` does

Inputs

- a directory of one-genome-per-file FASTA files
- a similarity output folder containing `similarities.tsv`
  -  Steps
1. Read all local FASTA files.
2. Read `similarities.tsv`.
3. Build genus and species clusters using the `70%` and `95%` thresholds - or set values 
4. Picks one representative per species cluster - using longest representative 
5. Write summary tables:
   - `cluster_assignments.tsv`
   - `genus_clusters.tsv`
   - `species_clusters.tsv`
   - `representatives.tsv`
6. Builds the custom DB assets:
   - `representatives.fasta`
   - `Bacteriophage_genomes.fasta`
   - `Bacteriophage_genomes.fasta.gz`
   - BLAST database files
   - `ICTV.msh`
   - `M.pa`- issue with this still 
   - `Mpa_pairwise_similarity.tsv`- ignore
   - `VMR.xlsx`

##### What `update_custom_taxmyphage_db.py` does

 Inputs

- an existing custom DB directory
- one or more new FASTA files or directories of FASTA files

 Steps

1. Load all genomes already known to the DB.
   - If `all_genomes.fasta` exists, use it.
   - Otherwise fall back to `representatives.fasta`.
2. Load the existing metadata from:
   - `all_genomes_metadata.tsv` if present
   - otherwise reconstruct enough metadata from `representatives.tsv` and `VMR.xlsx`
3. Combine old genomes and new genomes into one temporary FASTA.
4. Run `taxmyphage similarity` on the combined set.
5. Recomputes genus and species clusters from the new combined `similarities.tsv`.
6. Decide, for each new genome, whether it is:
   - `existing_species`
   - `new_species`
   - `new_genus`
7. Re-pick representatives from the updated species clusters.
8. Rebuild all DB assets from the updated representative set.
9. Write update outputs including:
   - `new_genome_report.tsv`
   - `all_genomes_metadata.tsv`
   - `last_update_similarities.tsv`

Naming behaviour

The updater preserves existing custom names where possible.

If a cluster contains no existing named genus/species, it allocates:

- `Custom_genus_XXX`
- `Custom_genus_XXX species_YYY`

If old custom names conflict because a previous database version had a bad split, the updater keeps the lowest-numbered compatible custom name as the canonical one.

Output files in the custom DB

- `all_genomes.fasta`
  All genomes currently represented in the maintained dataset.

- `all_genomes_metadata.tsv`
  Cluster assignment and naming metadata for every genome in the maintained dataset.

- `representatives.fasta`
  One representative genome per species cluster.

- `representatives.tsv`
  Representative metadata table.

- `cluster_assignments.tsv`
  Per-genome cluster assignments.

- `genus_clusters.tsv`
  Genus cluster membership.

- `species_clusters.tsv`
  Species cluster membership.

- `Bacteriophage_genomes.fasta` and `Bacteriophage_genomes.fasta.gz`
  FASTA used as the custom reference set for `taxmyphage`.

- `ICTV.msh`
  MASH sketch built from representative genomes staged in genus-specific subdirectories.

- `M.pa`
  Cached ordered-pair identity table among representatives - problematic and experimental 

- `Mpa_pairwise_similarity.tsv`
  Expanded pairwise table showing `idAB`, `idBA`, genome lengths, `simAB`, and `distAB`.

- `VMR.xlsx`
  Fake VMR workbook matching the representative set and custom names.

- `new_genome_report.tsv`
  Written by the updater to summarize how each newly added genome was classified.

## Typical usage

### Build an initial custom DB

```bash
python build_custom_taxmyphage_db.py --threads 4
```

### Update an existing custom DB with one genome

```bash
python update_custom_taxmyphage_db.py --db-dir custom_taxmyphage_db --input KCP_041.fa --threads 4
```

### Update with multiple new genomes

```bash
python update_custom_taxmyphage_db.py --db-dir custom_taxmyphage_db --input new_fastas_dir --threads 4
```

### Write the updated DB to a fresh directory instead of in place

```bash
python update_custom_taxmyphage_db.py \
  --db-dir custom_taxmyphage_db \
  --input new_fastas_dir \
  --threads 4 \
  --output-db-dir custom_taxmyphage_db_next
```

## 
