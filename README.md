# taxmyPHAGE

----------

Workflow to assign Taxonomy to a bacteriophage at the genus and species level.Briefly, the workflow will identify the most similar genomes in the set of currently classified ICTV genomes that are present in the VMR. 
Read about the VMR [here](https://ictv.global/vmr). It will compare the query genome against these genomes and run [VIRIDIC](https://doi.org/10.3390/v12111268) on the closest relatives. Interpret the output of VIRIDIC to determine if the phage falls within a current genus and or species. 



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


**VIRIDIC **

Requires viridic can be downloaded from here: http://rhea.icbm.uni-oldenburg.de/VIRIDIC/  


**MASH index**

A prebuilt MASH index of ICTV genomes. Can be downloaded from here https://millardlab-inphared.s3.climb.ac.uk/ICTV.msh

	wget  https://millardlab-inphared.s3.climb.ac.uk/ICTV.msh


**A database of genomes currently classified by the ICTV**

Can be created manually or download here  Bacteriophage Genomes https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz

	wget https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz
	gunzip Bacteriophage_genomes.fasta.gz

Create a blast database of these with 

	makeblastdb -in Bacteriophage_genomes.fasta -parse_seqids -dbtype nucl 


A copy of the VMR.xlsx - included here 


Run with 

	python tax_myPHAGE.py  --input UP30.fa

