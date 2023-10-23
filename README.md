# taxmyPHAGE

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/taxmyphage)](https://pypi.org/project/taxmyphage/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/pypi/v/taxmyphage)](https://pypi.org/project/taxmyphage/)
[![Conda](https://img.shields.io/conda/vn/bioconda/taxmyphage?style=plastic)](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/taxmyphage)
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


### A web version will be available soon. 

------

#### QUICK start and test

```
mamba create -n taxmyphage -c conda-forge -c bioconda taxymyphage
mamba activate taxmyphage

taxymyphage -i test.fna -t 4
```

if you are on macosx with a M1/M2 chip, you will need to install the following packages first

```
CONDA_SUBDIR=osx-64 mamba create -n taxmyphage -c conda-forge -c bioconda taxymyphage
mamba activate taxmyphage

taxymyphage -i test.fna -t 4
```

This should check the required software is installed and give a warning if not. It will also download the required fasta database and MASH file for comparison. These will be installed in the cloned tax_myPHAGE directory. If you download manually then please move them into tax_myPHAGE  directory.


Output of the test should have the following lines at the bottom 

 ![example](/img/example_result1.jpg)


----------

## Requirements 

It can be run on a standard laptop in a reasonable time. 


| **Input Files**                  | **Description**                                                                                                           | **Download Instructions**                                                |
|---------------------------------|---------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
| `ICTV.msh` (MASH index)         | Prebuilt MASH index of ICTV genomes. Download [here](https://millardlab-inphared.s3.climb.ac.uk/ICTV.msh) | Download: `wget https://millardlab-inphared.s3.climb.ac.uk/ICTV_2023.msh`  |
| `Bacteriophage_genomes.fasta.gz` | Database of ICTV-classified genomes. Download [here](https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz).                | - Download: `wget https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz`<br /> - Unzip: `gunzip Bacteriophage_genomes.fasta.gz`<br /> - Create blast db: `makeblastdb -in Bacteriophage_genomes.fasta -parse_seqids -dbtype nucl` |
| `VMR.xlsx`                      | Virus Metadata Resource (VMR) included.                                                                                             | - Download: `wget https://ictv.global/vmr/current`                                 |


------

### Install 

#### Conda 

```
mamba install -c conda-forge -c bioconda taxmyphage
```

if you are on macosx with a M1/M2 chip, you will need to install the following packages first

```
CONDA_SUBDIR=osx-64 mamba install -c conda-forge -c bioconda taxmyphage
```

Bioconda doesn't support osx-arm64 yet.


#### pip

```
pip install taxmyphage
```

If you installing by pip, don't forget to install a working version of [mash](https://mash.readthedocs.io/en/latest/) and [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


### Run with 

#### Quick start

```
taxmyphage -i input.fasta
```

#### Options 

```
optional arguments:
  -h, --help            show this help message and exit
  -V, --version         Show the version number and exit.
  -v, --verbose

General options:
  -i [FASTA_FILE ...] [[FASTA_FILE ...] ...], --input [FASTA_FILE ...] [[FASTA_FILE ...] ...]
                        Path to an input fasta file(s), or directory containing fasta files. The fasta file(s) could contain multiple phage genomes. They can be compressed (gzip). If a directory is given the
                        expected fasta extentions are ['fasta', 'fna', 'fsa', 'fa'] but can be gzipped. (Required)
  -o OUTPUT, --output OUTPUT
                        Path to the output directory. (Default is current directory)
  -p PREFIX, --prefix PREFIX
                        will add the prefix to results and summary files that will store results of MASH and comparision to the VMR Data produced by ICTV combines both sets of this data into a single csv file.
                        Use this flag if you want to run multiple times and keep the results files without manual renaming of files. (Default no prefix)
  -d DIST, --distance DIST
                        Will change the mash distance for the intial seraching for close relatives. We suggesting keeping at 0.2 If this results in the phage not being classified, then increasing to 0.3 might
                        result in an output that shows the phage is a new genus. We have found increasing above 0.2 does not place the query in any current genus, only provides the output files to demonstrate it
                        falls outside of current genera. (Default is 0.2)
  --no-figures          Use this option if you don't want to generate Figures. This will speed up the time it takes to run the script - but you get no Figures. (By default, Figures are generated)
  -t THREADS, --threads THREADS
                        Maximum number of threads that will be used by BLASTn. (Default is 1)

Options related to the database:
  --db_folder FOLDER_PATH
                        Path to the database directory where the databases are stored. (Default is /Users/rdenise/Documents/Scripts/tax_myPHAGE_github_version/taxmyphage/database)

Install folder and database options:
  --install             Use this option to download and install the databases. (Default is False)

Executable options, if not in PATH:
  --blastcmd BLASTDBCMD   Path to the blastn executable (default: blastdbcmd)
  --blastn BLASTN       Path to the blastn executable (default: blastn)
  --makeblastdb MAKEBLASTDB
                        Path to the blastn executable (default: makeblastdb)
  --mashexe MASHEXE     Path to the MASH executable (default: mash)
```
----------



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



#### Output files 

```
[output_folder]                                          <- General output folder
└── [genome query_id]                                    <- Results output for the genome query_id
    ├── Output_of_taxonomy.csv                           <- Output of the taxonomy classification
    ├── Summary_file.txt                                 <- Summary of the analysis (summarises what was printed to screen )
    ├── heatmap.pdf                                      <- Heatmap of the similarity to the closest relatives (pdf)
    ├── heatmap.png                                      <- Heatmap of the similarity to the closest relatives    (png)
    ├── heatmap.svg                                      <- Heatmap of the similarity to the closest relatives (svg)
    ├── known_taxa.fa                                    <- Fasta file of the closest relatives
    ├── mash.txt                                         <- Output of the MASH analysis
    ├── query.fasta                                      <- Input fasta file
    ├── similarities.tsv                                 <- Similarities to the closest relatives
    ├── top_right_matrix.tsv                             <- Top right matrix of similarity to closest relatives (same as heatmap)
    └── viridic                                          <- VIRIDIC-like analysis
        ├── viridic_in.fa                                <- Input fasta file
        ├── viridic_in.fa.blastn_vs2_self.tab.gz         <- Blastn output of the input fasta file against itself
        ├── viridic_in.fa.genus_species_clusters.tsv     <- Clusters of the closest relatives
        ├── viridic_in.fa.ndb                            <- Blastn database of the closest relatives
        ├── viridic_in.fa.nhr
        ├── viridic_in.fa.nin
        ├── viridic_in.fa.njs
        ├── viridic_in.fa.not
        ├── viridic_in.fa.nsq
        ├── viridic_in.fa.ntf
        └── viridic_in.fa.nto
```

##### Output files explained

- **Summary_file.txt** - Summarises what was printed to screen 

```
Query sequence header was:test1 
    
    
Query sequence can be classified within a current genus and represents a new species, it is in:
    
Class:Caudoviricetes    Family: Not Defined Yet    Subfamily:Vequintavirinae    Genus:Certrevirus    Species:name_your_species
```

- **Output_of_taxonomy.csv** - Provides Cluster and Species numbers for you query phage, merged with data from the VMR for the closest relatives to you query

- ***pdf, *svg, *png**  - image files of top right matrix of similarity to closest currently classified phages 

 ![HeatMap](/img/heatmap.jpg)
  
    
