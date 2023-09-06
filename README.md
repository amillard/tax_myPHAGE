# taxmyPHAGE

----------

Workflow to assign Taxonomy to a bacteriophage at the genus and species level.Briefly, the workflow will identify the most similar genomes in the set of currently classified ICTV genomes that are present in the VMR. 
Read about the VMR [here](https://ictv.global/vmr). It will compare the query genome against these genomes and run a [VIRIDIC](https://doi.org/10.3390/v12111268)-**like analysis** on the closest relatives. Interpret the output of VIRIDIC-like analysis to determine if the phage falls within a current genus and or species. 



Designed for:

- Individual complete phage genomes 


What it wont do:
 
- Work on multiple  genomes in a file 
- Metagenomic samples 
- Eukaryotic viruses
- RNA phages - it will give a result - not necessarily the correct one 
- ssDNA phages - again a result but likely not accurate 
- Classify a phage into a new family 




## Requirements 

----------

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

	python tax_myPHAGE.py  --input UP30.fa 

or 

	python taxmyPHAGE.py ---h to get full options 



	-h, --help            show this help message and exit


	  -v, --verbose
	  -t THREADS, --threads THREADS
	                        Maximum number of threads that will be used

	  -i IN_FASTA, --input IN_FASTA
	                        Path to an input fasta file

	 
	  -d DIST, --distance DIST
	                        Will change the mash distance for the intial seraching for close relatives. We suggesting keeping
	                        at 0.2 If this results in the phage not being classified, then increasing to 0.3 might result in
	                        an output that shows the phage is a new genus. We have found increasing above 0.2 does not place
	                        the query in any current genus, only provides the output files to demonstrate it falls outside of
	                        current genera

	  --Figures {T,F}       Specify 'T' or 'F' to produce Figures. Using F will speed up the time it takes to run the script -
	                        but you get no Figures. Default is with Figures 


----------



#### Basic benchmarks 




| Genus | Number of genomes in Genera|Time(H:M:S)
| ------------- | ------------- |-------
|Cheoctovirus |96|00:07:44
|Tequatrovirus|83|00:26:19|
|Peduovirus |27|00:00:23|
|Warwickvirus|18|00:00:18|
|Pseudotevenvirus|9|0:01:15|
|Changmaivirus|2|0:00:17
|Stompvirus|1|0:00:16



