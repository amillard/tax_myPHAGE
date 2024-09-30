# taxmyPHAGE

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/taxmyphage)](https://pypi.org/project/taxmyphage/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/pypi/v/taxmyphage)](https://pypi.org/project/taxmyphage/)
[![Conda](https://img.shields.io/conda/vn/bioconda/taxmyphage)](https://anaconda.org/bioconda/taxmyphage)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/taxmyphage)](https://pepy.tech/project/taxmyphage)

----------

Developed by [Andrew Millard](https://github.com/amillard) , [Thomas Sicheritz-Ponten](https://github.com/tsp-kucbd) and [Remi Denise](https://github.com/rdenise)

Script to assign taxonomy to a bacteriophage at the genus and species level. It will identify the most similar genomes in the set of currently classified ICTV genomes that are present in the VMR. 
Read about the VMR [here](https://ictv.global/vmr). It will compare the query genome against these genomes and run a [VIRIDIC](https://doi.org/10.3390/v12111268)-**like analysis** on the closest relatives. Interpret the output of VIRIDIC-like analysis to determine if the phage falls within a current genus and or species. It does not run VIRIDIC, but utilises the same formula for comparison of genomes.  The input is a single genome sequence. The remainder of the analysis is automated 

----------

Designed for:

- Individual complete phage genomes 

What it will do:

- Classify a dsDNA phage genomes at the Genus and or species level against ICTV genomes 
- Tell you if your genome represents a new genus 
- Use current ICTV cutoffs for Genera and Species 
- accept multiple inputs at the same time

----------

What it wont do:
 
- Metagenomic samples 
- Eukaryotic viruses
- RNA phages - it will give a result - not necessarily the correct one 
- ssDNA phages - again a result but likely not accurate 
- Classify a phage into a new family 
- Compare against every single phage genome in Genbank. It is designed for classification , so compares against currently classified phages.


*Webversion is [here](https://ptax.ku.dk/)*

------


## QUICK start and test

```
mamba create -n taxmyphage -c conda-forge -c bioconda taxmyphage
mamba activate taxmyphage

# If databases not installed, install them
taxmyphage install

# Run taxmyphage
taxmyphage run -i test.fna -t 4
```

if you are on macosx with a M1/M2 chip, you will need to install the following packages first

```
CONDA_SUBDIR=osx-64 mamba create -n taxmyphage -c conda-forge -c bioconda taxmyphage
mamba activate taxmyphage

# If databases not installed, install them
taxmyphage install

# Run taxmyphage
taxmyphage run -i test.fna -t 4
```

This should check the required software is installed and give a warning if not. It will also download the required fasta database and MASH file for comparison. These will be installed in the cloned tax_myPHAGE directory. If you download manually then please move them into tax_myPHAGE  directory.


Output of the test should have the following lines at the bottom 

<img src="/img/example_result.png" width="70%" height="70%">


----------


# Manual update for new VMR 

With the relase of a new VMR the databases for taxmyphage need to be updated if you installed prior to the 17th May 2024. For new installs you will also need to do the following to use the current [VMR](https://ictv.global/vmr).  While we work on an automated solution, the following can be done. 

Install taxmyphage using any of the methods listed below. 

Then update the database files with the links below.  

New Genome database [here](https://millardlab-taxmyphage.s3.climb.ac.uk/Bacteriophage_genomes_MSL39_v1.fasta.gz)
New MASH database [here](https://millardlab-taxmyphage.s3.climb.ac.uk/ICTV_MSL39v1.msh) 

Find the files Bacteriophages_genomes.fasta.gz in the conda installation  eg 
```
find ~/miniconda/ -name "Bacteriophage_genomes*"
```

Do the same for ICTV.msh

```
> find ~/miniconda*/ -name "ICTV*"
```

and the VMR.xlsx file 

```
find ~/miniconda*/ -name "VMR*"
```

Delete these files and then replace with files downloaded above and rename them.

my files are in

/home/andrew/miniconda3/envs/taxmyphage/lib/python3.12/site-packages/taxmyphage/database/ICTV_2023.msh

```
mv Bacteriophage_genomes_MSL39_v1.fasta.gz  Bacteriophage_genomes.fasta.gz
```

```
mv ICTV_MSL39v1.msh ICTV.msh
```

Or just create softlinks to the oringinal name.  

It will now run with the latest version of the ICTV taxonomy 


### Change to VMR format September 29th 2024

An update to the format of VMR file produced by the ICTV and released ~20 th  September 2024 , has different format to their previous files. If you install after this date, it will download the latest version of the VMR (VMR_MSL39_v2.xlsx) and cause taxMyPhage to not run correctly. We are working on a solution to correctly format the newest VMR format to our input needs for taxMyPhage. In the mean time the easiest solution is to manually download and install VMR_MSL39_v1.xlsx rather than VMR_MSL39_v2.xlsx and rename the file accordingly , as above. 


## Requirements 

It can be run on a standard laptop in a reasonable time. 


| **Input Files**                  | **Description**                                                                                                           | **Download Instructions**                                                |
|---------------------------------|---------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
| `ICTV.msh` (MASH index)         | Prebuilt MASH index of ICTV genomes. Download [here](https://millardlab-inphared.s3.climb.ac.uk/ICTV.msh) | Download: `wget https://millardlab-inphared.s3.climb.ac.uk/ICTV_2023.msh`  |
| `Bacteriophage_genomes.fasta.gz` | Database of ICTV-classified genomes. Download [here](https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz).                | - Download: `wget https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz`<br /> - Unzip: `gunzip Bacteriophage_genomes.fasta.gz`<br /> - Create blast db: `makeblastdb -in Bacteriophage_genomes.fasta -parse_seqids -dbtype nucl` |
| `VMR.xlsx`                      | Virus Metadata Resource (VMR). Download [here](https://ictv.global/vmr/current).                                                                                             | - Download: `wget https://ictv.global/vmr/current`                                 |


------

## Install 

> [!IMPORTANT]  
> tax_myPHAGE requires [MASH](https://mash.readthedocs.io/en/latest/) >=2.3 and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) >=2.14.0.  
> You need to install mash and NCBI BLAST by yourself (except if you install macsyfinder via *conda/mamba*).  
> The other dependencies are managed by the python package manager *pip*.  

### Conda 

```
mamba install -c conda-forge -c bioconda taxmyphage
```

if you are on macosx with a M1/M2 chip, you will need to install the following packages first

```
CONDA_SUBDIR=osx-64 mamba install -c conda-forge -c bioconda taxmyphage
```

Bioconda doesn't support osx-arm64 yet.


### Pypi

```
pip install taxmyphage
```

If you installing by pip, don't forget to install a working version of [MASH](https://mash.readthedocs.io/en/latest/) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

### Source

Alternatively, a development version of `taxmyphage` can be install from github.

```
git clone https://github.com/amillard/tax_myPHAGE
```

either you can install the taxmyphage and its dependencies via pip

```
cd tax_myPHAGE
pip install -e .
taxmyphage -h
```

or you install following dependencies yourself

```
biopython >= 1.81
pandas >= 2.1.1
seaborn >= 0.13
wget >= 3.2
scipy >= 1.11.3
tqdm >= 4.66.1
openpyxl >= 3.1.2
networkx >= 3.1
icecream >= 2.1.3
```

and run it via

```
pip install -e --no-dependencies .
taxmyphage -h
```

or

```
tax_myPHAGE/taxmyphage/bin/taxmyphage.py -h
```


## Modules

### Install

Allows to install the databases before running the `Run` module.

```
taxmyphage install
```

#### Options 

```
usage: taxmyphage install [-h] [-v] [-V] [-db FOLDER_PATH] [--makeblastdb MAKEBLASTDB]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Show verbose output. (For debugging purposes)
  -V, --version         Show the version number and exit.

Databases options:
  -db FOLDER_PATH, --db_folder FOLDER_PATH
                        Path to the database directory where the databases are stored. (Default is /Users/rdenise/Documents/Scripts/tax_myPHAGE/taxmyphage/database)

Install options:
  --makeblastdb MAKEBLASTDB
                        Path to the blastn executable (default: makeblastdb)
```
----------
### Run

```
taxmyphage run -i input.fasta
```

#### Options 

```
usage: taxmyphage run [-h] -i [FASTA_FILE ...] [[FASTA_FILE ...] ...] [-o OUTPUT] [-p PREFIX] [-t THREADS] [-d DIST] [--mash MASH] [--blastdbcmd BLASTDBCMD] [--blastn BLASTN] [--makeblastdb MAKEBLASTDB]
                      [--no-figures] [-v] [-V] [-db FOLDER_PATH]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Show verbose output. (For debugging purposes)
  -V, --version         Show the version number and exit.

General options:
  -i [FASTA_FILE ...] [[FASTA_FILE ...] ...], --input [FASTA_FILE ...] [[FASTA_FILE ...] ...]
                        Path to an input fasta file(s), or directory containing fasta files. The fasta file(s) could contain multiple phage genomes. They can be compressed (gzip). If a directory is given the
                        expected fasta extentions are ['fasta', 'fna', 'fsa', 'fa'] but can be gzipped. (Required)
  -o OUTPUT, --output OUTPUT
                        Path to the output directory. (Default is current directory)
  -p PREFIX, --prefix PREFIX
                        Will add the prefix to results and summary files that will store results of MASH and comparision to the VMR Data produced by ICTV combines both sets of this data into a single csv file.
                        Use this flag if you want to run multiple times and keep the results files without manual renaming of files. (Default no prefix)
  -t THREADS, --threads THREADS
                        Maximum number of threads that will be used by BLASTn. (Default is 1)

MASH options:
  -d DIST, --distance DIST
                        Will change the mash distance for the intial seraching for close relatives. We suggesting keeping at 0.2 If this results in the phage not being classified, then increasing to 0.3 might
                        result in an output that shows the phage is a new genus. We have found increasing above 0.2 does not place the query in any current genus, only provides the output files to demonstrate it
                        falls outside of current genera. (Default is 0.2)
  --mash MASH           Path to the MASH executable (default: mash)
  --blastdbcmd BLASTDBCMD
                        Path to the blastn executable (default: blastdbcmd)

Similarity options:
  --blastn BLASTN       Path to the blastn executable (default: blastn)
  --makeblastdb MAKEBLASTDB
                        Path to the blastn executable (default: makeblastdb)
  --no-figures          Use this option if you don't want to generate Figures. This will speed up the time it takes to run the script - but you get no Figures. (By default, Figures are generated)

Databases options:
  -db FOLDER_PATH, --db_folder FOLDER_PATH
                        Path to the database directory where the databases are stored. (Default is /Users/rdenise/Documents/Scripts/tax_myPHAGE/taxmyphage/database)
```


#### Indicative run time  

The time to classify a phage will depend on the number of hits and number of phages currently classified within a particular genus. The more species within a genus, the longer the time for classification. The numbers below are from running on a 16 core server. We have been running the process on a MAC book and Windows laptop in reasonable time periods. 



| Genus | Number of genomes in Genera|Time(H:M:S)
| ------------- | ------------- |-------
|Cheoctovirus |96|00:07:44
|Tequatrovirus|83|00:26:19|
|Peduovirus |27|00:00:23|
|Warwickvirus|18|00:00:18|
|Pseudotevenvirus|9|0:01:15|
|Changmaivirus|2|0:00:17
|Stompvirus|1|0:00:16



#### Output files for the `Run` module

```
[output_folder]                                          <- General output folder
├── Summary_taxonomy.tsv                                 <- Summary of the analysis for all the genomes (summarises what was printed to screen)
└── Results_per_genome                                   <- Folder containing the results for each genome
  └── [genome query_id]                                  <- Results output for the genome query_id
      ├── Output_of_taxonomy.tsv                         <- Output of the taxonomy classification
      ├── Summary_file.txt                               <- Summary of the analysis (summarises what was printed to screen)
      ├── heatmap.pdf                                    <- Heatmap of the similarity to the closest relatives (pdf)
      ├── heatmap.png                                    <- Heatmap of the similarity to the closest relatives    (png)
      ├── heatmap.svg                                    <- Heatmap of the similarity to the closest relatives (svg)
      ├── known_taxa.fa                                  <- Fasta file of the closest relatives
      ├── mash.txt                                       <- Output of the MASH analysis
      ├── query.fasta                                    <- Input fasta file
      ├── similarities.tsv                               <- Similarities to the closest relatives
      ├── top_right_matrix.tsv                           <- Top right matrix of similarity to closest relatives (same as heatmap)
      └── pmv                                            <- Clustering on genomic similarity analysis
          ├── pmv_in.fa                                  <- Input fasta file
          ├── pmv_in.fa.blastn_vs2_self.tab.gz           <- Blastn output of the input fasta file against itself
          ├── pmv_in.fa.genus_species_clusters.tsv       <- Clusters of the closest relatives
          ├── pmv_in.fa.ndb                              <- Blastn database of the closest relatives
          ├── pmv_in.fa.nhr
          ├── pmv_in.fa.nin
          ├── pmv_in.fa.njs
          ├── pmv_in.fa.not
          ├── pmv_in.fa.nsq
          ├── pmv_in.fa.ntf
          └── pmv_in.fa.nto
```

##### Output files explained

- **Summary_taxonomy.tsv** - Summarises what was printed to screen for all the genomes

| Query sequence header | Realm | Kingdom | Phylum | Class | Order | Family | Subfamily | Genus | Species | Full classification | Message |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| MZ130489 | Duplodnaviria | Heunggongvirae | Uroviricota | Caudoviricetes | Crassvirales | Crevaviridae | Coarsevirinae | Junduvirus | Junduvirus communis | r__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes;o__Crassvirales;f__Crevaviridae;sf__Coarsevirinae;g__Junduvirus;s__Junduvirus communis | Current ICTV taxonomy and clustering on genomic similarity algorithm output appear to be consistent at the genus level |
| newGenus_phage | Unknown | Unknown | Unknown | Unknown | Unknown | Unknown | Unknown | New_genus | New_genus new_species | r__Unknown;k__Unknown;p__Unknown;c__Unknown;o__Unknown;f__Unknown;sf_Unknown;g__New_genus;s__New_genus new_species | Query is a new genus and species. You could try running again with if you larger distance |
| test1 | Duplodnaviria | Heunggongvirae | Uroviricota | Caudoviricetes | Not Defined Yet | Not Defined Yet | Vequintavirinae | Certrevirus | Certrevirus name_your_species | r__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes;o__;f__;g__Certrevirus;s__Certrevirus name_your_species | The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation |

- **Summary_file.txt** - Summarises what was printed to screen 

```
Query sequence header was:test1 
    
    
Query sequence can be classified within a current genus and represents a new species, it is in:
    
Class:Caudoviricetes    Family: Not Defined Yet    Subfamily:Vequintavirinae    Genus:Certrevirus    Species:name_your_species
```

- **Output_of_taxonomy.csv** - Provides Cluster and Species numbers for you query phage, merged with data from the VMR for the closest relatives to you query

- ***.pdf, *.svg, *.png**  - image files of top right matrix of similarity to closest currently classified phages 

<img src="/img/heatmap.png" width="60%" height="60%">

----------
### Similarity

```
taxmyphage similarity -i input.fasta
```

#### Options 

```
usage: taxmyphage similarity [-h] -i [FASTA_FILE ...] [[FASTA_FILE ...] ...] [-o OUTPUT] [-p PREFIX] [-t THREADS] [--reference REFERENCE] [--blastn BLASTN] [--makeblastdb MAKEBLASTDB] [--no-figures] [-v] [-V]
                          [-db FOLDER_PATH]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Show verbose output. (For debugging purposes)
  -V, --version         Show the version number and exit.

General options:
  -i [FASTA_FILE ...] [[FASTA_FILE ...] ...], --input [FASTA_FILE ...] [[FASTA_FILE ...] ...]
                        Path to an input fasta file(s), or directory containing fasta files. The fasta file(s) could contain multiple phage genomes. They can be compressed (gzip). If a directory is given the
                        expected fasta extentions are ['fasta', 'fna', 'fsa', 'fa'] but can be gzipped. (Required)
  -o OUTPUT, --output OUTPUT
                        Path to the output directory. (Default is current directory)
  -p PREFIX, --prefix PREFIX
                        Will add the prefix to results and summary files that will store results of MASH and comparision to the VMR Data produced by ICTV combines both sets of this data into a single csv file.
                        Use this flag if you want to run multiple times and keep the results files without manual renaming of files. (Default no prefix)
  -t THREADS, --threads THREADS
                        Maximum number of threads that will be used by BLASTn. (Default is 1)

Comparison options:
  --reference REFERENCE
                        Path to the reference database file. Input file will be used as query against it. If not provided, input will be compare against itself. If you use reference no figure is generated.
                        (Default is '')

Similarity options:
  --blastn BLASTN       Path to the blastn executable (default: blastn)
  --makeblastdb MAKEBLASTDB
                        Path to the blastn executable (default: makeblastdb)
  --no-figures          Use this option if you don't want to generate Figures. This will speed up the time it takes to run the script - but you get no Figures. (By default, Figures are generated)

Databases options:
  -db FOLDER_PATH, --db_folder FOLDER_PATH
                        Path to the database directory where the databases are stored. (Default is /Users/rdenise/Documents/Scripts/tax_myPHAGE/taxmyphage/database)
```

#### Output files for the `similarity` module

```
[output_folder]                          <- General output folder
├── heatmap.pdf                          <- Heatmap of the similarity to the closest relatives (pdf)
├── heatmap.png                          <- Heatmap of the similarity to the closest relatives    (png)
├── heatmap.svg                          <- Heatmap of the similarity to the closest relatives (svg)
├── pmv.fasta                            <- Input fasta file
├── pmv.fasta.blastn_vs2_self.tab.gz     <- Blastn output of the input fasta file against itself
├── pmv.fasta.genus_species_clusters.tsv <- Clusters of the closest relatives
├── pmv.fasta.ndb                        <-  Blastn database of the closest relatives
├── pmv.fasta.nhr
├── pmv.fasta.nin
├── pmv.fasta.njs
├── pmv.fasta.not
├── pmv.fasta.nsq
├── pmv.fasta.ntf
├── pmv.fasta.nto
├── similarities.tsv                     <- Similarities to the closest relatives
└── top_right_matrix.tsv                 <- Top right matrix of similarity to closest relatives (same as heatmap)
```

##### Output files explained

- ***.pdf, *.svg, *.png**  - image files of top right matrix of similarity to closest currently classified phages 

<img src="/img/heatmap.png" width="60%" height="60%">
