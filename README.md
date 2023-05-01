# taxamyphage
Workflow to assign Taxonomy to a bacteriophage 

###VIRIDIC
Requires viridic can be downloaded from here: http://rhea.icbm.uni-oldenburg.de/VIRIDIC/  

Also requires MASH index of ICTV genomees. Can be downloaded from here https://millardlab-inphared.s3.climb.ac.uk/ICTV.msh

	`wget  https://millardlab-inphared.s3.climb.ac.uk/ICTV.msh  `


Also requires the ICTV classified Bacteriophage Genomes https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz

	`wget https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz  `
	`gunzip  https://millardlab-inphared.s3.climb.ac.uk/Bacteriophage_genomes.fasta.gz`

Create a blast database of these with 

	` makeblastdb -in Bacteriophage_genomes.fasta -parse_seqids -dbtype nucl` 


Run with 

	`python python tax_myPHAGE.py  --input UP30.fa`

