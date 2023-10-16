# taxmyPHAGE

----------

Developed by A Millard , [Thomas Sicheritz-Ponten](https://github.com/tsp-kucbd) and [Remi Denise](https://github.com/rdenise)

Script to assign taxonomy to a bacteriophage at the genus and species level. It will identify the most similar genomes in the set of currently classified ICTV genomes that are present in the VMR. 
Read about the VMR [here](https://ictv.global/vmr). It will compare the query genome against these genomes and run a [VIRIDIC](https://doi.org/10.3390/v12111268)-**like analysis** on the closest relatives. Interpret the output of VIRIDIC-like analysis to determine if the phage falls within a current genus and or species. It does not run VIRIDIC, but utilises the same formula for comparison of genomes.  The input is a single genome sequence. The remainder of the analysis is automated 



Designed for:

- Individual complete phage genomes 

What it will do:

- Classify a dsDNA phage genomes at the Genus and or species level against ICTV genomes 
- Tell you if your genome represents a new genus 
- Use current ICTV cutoffs for Genera and Species 
- accept multiple inputs at the same time (not a multifasta)



What it wont do:
 
- Work on multiple  genomes in a file 
- Metagenomic samples 
- Eukaryotic viruses
- RNA phages - it will give a result - not necessarily the correct one 
- ssDNA phages - again a result but likely not accurate 
- Classify a phage into a new family 
- Compare against every single phage genome in Genbank. It is designed for classification , so compares against currently classified phages.


- ### A web version will be available soon. 

------

#### QUICK start and test


	git clone https://github.com/amillard/tax_myPHAGE

	cd tax_myPHAGE

	mamba install  -c conda-forge -c bioconda biopython pandas icecream networkx tqdm openpyxl matplotlib scipy

	python tax_myPHAGE.py -i test.fna -t 8 


This should check the required software is installed and give a warning if not. It will also download the required fasta database and MASH file for comparison. These will be installed in the cloned tax_myPHAGE directory. If you download manually then please move them into tax_myPHAGE  directory.


Output of the test should have the following lines at the bottom 

 ![example](/img/example_result1.jpg)





## Requirements 

----------

It can be run on a standard laptop in a reasonable time. 




### MASH  

A working version of [mash](https://mash.readthedocs.io/en/latest/) for install instructions


**MASH index**

A prebuilt MASH index of ICTV genomes. Can be downloaded from here https://millardlab-inphared.s3.climb.ac.uk/ICTV.msh

	wget  https://millardlab-inphared.s3.climb.ac.uk/ICTV_2023.msh

Will attempt to install automatically if you haven't downloaded in advance of running 


### **A database of genomes currently classified by the ICTV**

Can be created manually or download here [Bacteriophage Genomes](https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz)

	wget https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz
	gunzip Bacteriophage_genomes.fasta.gz

Should be included in the cloned tax_myPHAGE directory 



Create a blast database of these with 

	makeblastdb -in Bacteriophage_genomes.fasta -parse_seqids -dbtype nucl 


Again it will attempt to download and install these for you if they havent been installed in advance 

----------

### VMR


A copy of the VMR.xlsx - included here 

Again will download a version if none is not detected 



------

### Install python modules 


	mamba install biopython pandas icecream networkx tqdm openpyxl matplotlib





### Run with 

	python tax_myPHAGE.py  --input input.fa 

or 

	python taxmyPHAGE.py ---h to get full options 



	-h, --help            show this help message and exit


	  -v, --verbose
	  -t THREADS, --threads THREADS
	                        Maximum number of threads that will be used

	  -i IN_FASTA, --input IN_FASTA
	                        Path to an input fasta file or multiple input files eg *fasta or genome1.fna genome2.fa  
							or a directory with fasta files in -i my_genome_dir 
						    or mix of a directory and files -i my_genome_dir genome1.fna
							Will recognise *fasta *fna *fsa *fa and gz of these files 
							DO NOT USE -p option when specifying multiple files. All output will have that prefix!
							

	 
	  -d DIST, --distance DIST
	                        Will change the mash distance for the intial seraching for close relatives. We suggesting keeping
	                        at 0.2 If this results in the phage not being classified, then increasing to 0.3 might result in
	                        an output that shows the phage is a new genus. We have found increasing above 0.2 does not place
	                        the query in any current genus, only provides the output files to demonstrate it falls outside of
	                        current genera

	  --Figures {T,F}       Specify 'T' or 'F' to produce Figures. Using F will speed up the time it takes to run the script -
	                        but you get no Figures. Default is with Figures 


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






##### Output files 


- **Summary_file.txt** - summarises what was printed to screen 


eg

Query sequence header was:test1 
	
	
Query sequence can be classified within a current genus and represents a new species, it is in:
	
Class:Caudoviricetes	Family: Not Defined Yet	Subfamily:Vequintavirinae	Genus:Certrevirus	Species:name_your_species

---------

- **Output_of_taxonomy.csv** - Provides Cluster and Species numbers for you query phage, merged with data from the VMR for the closest relatives to you query

- ***pdf, *svg, *jpg**  - image files of top right matrix of similarity to closest currently classified phages 



 ![HeatMap](/img/heatmap.jpg)
  
    
