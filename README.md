# NaDA
Nanopore Data Analysis(NaDA) is a bachelor internship program, focussing on analysing and visualization of ctDNA for a p53 mutant database. NaDA provides a range of python scripts able to analyse single nucleotide polymorphisms (SNPs) in both fastq and vcf files derived from nanopore sequencing.

# Documentation
further documentation about NaDA can be found [here](https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/index.html) and for more background information my literature studie reveals more about the importance of ctDNA analysis.

# Installation and Usage
The following scripts in Nada/Scripts can be used to analyse fastq and vcf files:
1) fastq_qualityscore_analyser.py: analyses fastq files for both scores by sequence and cluster
2) fastq_meanscore_plot.py: analyses save data from fastq_qualityscore_analyser.py to determine and plot run meanscore
3) vcf_snp_variant_analyser.py: analyses vcf files for occurence of variance in both specific location and a specific sequence
4) vcf_position_mutation.py: analyses csv files from vcf_snp_variant_analyser.py to determine percentage of specific mutation on each location. Can analyse mulitple csv files simultaneously. 

Scripts are highly customizable in terms of sequence size selection and backbone selection

Argument requirements can be identified using help statement. Scripts have been tested on python 3.6.1 and this version of python is therefore advised. Further package requirements can be found in Nada/requirements.txt

Project Nada can also be run on a High-Performance Computing system. The HPC facility consists of 1544 cores and 600TB of High-Performance storage. The HPC facility runs on CentOS Linux and provides a batch-wise queueing system with a few head nodes and many compute nodes for submitting and running many computational tasks in parallel. For our project HPC facility in the Utrecht Science Park were used. In order to use these scripts on the HPC, two run scripts have been written. run_NaDA.py is used to submit vcf_snp_variant_analyser.py to the HPC facility and run_fastq.py is used to submit fastq_qualityscore_analyser.py to the HPC facility. !MPORTANT! In order to use these run scripts, make sure python 3.6.1 is selected. Furthermore, cmd0 activates a virtual enviroment which should include earlier described python packages. the qsub_cmd command line can be customized to own preference. Currently, 24 hour run time, 32GB of ram and 8 threads have been selected as ell as my email. 
