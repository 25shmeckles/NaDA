.. Internship report documentation master file, created by
   sphinx-quickstart on Wed Dec  6 10:37:32 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Data analyzing of ctDNA from Nanopore Sequencing with Python
============================================================
**Cyclomics data analysis: identifying and analyzing alterations in ctDNA for data filtering**

*Student:			Douwe Spaanderman (5483352)*

*Supervisor: 			Alessio Marcozzi*

*Group Kloosterman		Medical Genetics - UMC Utrecht*


Introduction
------------
Tumor genotyping is essential for improving the clinical outcome for many cancer patients. Not only is it important for tailoring drugs specific targets for the genetic landscape of the tumor. It is also vital for better diagnostics as well as follow up monitoring. Furthermore, increasing knowledge of tumor tissue genotyping can enhance our understanding of the molecular mechanisms of cancer. Today, tumor genotyping is conducted using local biopsies of tumor tissue. This method of achieving the genetic landscape is however heavily limited mostly due to its invasive nature. Therefore, tissue biopsies can only be conducted a limited amount of times, enhancing difficulties with tumor heterogeneity and the fluctuation of the genetic landscape of the tumor in time(1). Altogether these limitations require a new none invasive method to be able to tackle these difficulties.

Just recently, circulating tumor DNA (ctDNA) has been suggested as a new innovating method for identifying the tumor genetic landscape. These DNA fragments are derived from necrotic tumor cells and subsequently further fragmented by macrophages(2). ctDNA makes up for a portion of the whole cell-free DNA (cfDNA), increasingly so upon heavy tumor burden. Theoretically, ctDNA provides a better alternative for tumor genotyping, because of its noninvasive nature, ability to be used unlimited times, it accounts for most of both intratumoral and intermetastatic heterogeneity and has been shown to identify multiple types of alterations. However, ctDNA encounters some difficulties, such as low allelic fractions of ctDNA in patient's blood(3), especially on low tumor burden. Also, discrimination between cfDNA and ctDNA has seemed to be problematic. These problems hinder researchers to develop high specific and sensitive clinical tests.

The Kloosterman group(UMC Utrecht) uses a new and innovating method called Cyclomics to identify alterations in ctDNA (Figure 1). Currently research upon this method is conducted for p53 mutation in patients with head and neck carcinoma's, but upon achieving a sensitive and specific test, Cyclomics could provide for identifying multiple driver mutations in ctDNA(4). Cyclomics uses rolling circle amplification to increase low allelic fractions of mutated ctDNA fragments. This library preparation for Nanopore sequencing increases its sensitivity. However, as with all sequencing methods, sequencing isn't flawless and covers real biological mutations but also asymmetric DNA errors, PCR and sequencing errors(5). Therefore, including multiple filtering steps upon the data achieved from the Nanopore sequencer is essential for identifying the real biological mutations and thus improving specificity and sensitivity of the clinical tests. In this report, data analyzing will be covered with Python, identifying regions with high and low mutational value as well as the higher occurrence of certain single nucleotide polymorphisms (SNPs).
 
.. figure::  C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Figure_workflow.png
   :scale:   70%
   :align:   center

   Figure 1: Workflow, NOTE: made in paint, this is the test version

**General workflow**

Results and discussion
----------------------
Fastq files
+++++++++++

#Where are both the Fastq and vcf files derived from. Are they from healthy individuals, how many persons, ect.

**Quality score distribution**. Fastq are text-based files derived from sequencers such as Illumina and Nanopore that contain a biological sequence and for every base a quality score.
this quality score represents the sequencers ability to call a base as true, in other words, the score says something about how sure the sequencer is that the right base is called. 
if the score for the base is low, it's most likely the right base, while if the score is high the machine is unsure if it's the real base. Analyzing quality scores on sequences
is important to be able to understand which regions of bases is difficult to get base called. More importantly, a mutation called in a region with low quality score is therefore more 
likely to be a real mutation in contrast to a mutation that occurs through wrong base calling in area's with high quality scores. Clearifying which regions are high score is important
for achieving tests with high specificity, as these regions will not be selected for variant calling. Likewise, regions with low score are good targets for identifying mutations 
and thus are able to increase tests sensitivity. In order to analyse fastq files the following script was written to strip the files:: 

	def parse_fasta_file_error(sequence_file):
		data = {}
			with open(sequence_file, 'r') as f:
        		id_ = False
       			sequence = False
       			for line in f:
       				if line.strip() in ['+','\n']:
               				continue
           			elif line.startswith('@'):
                			if not id_: #if id_ == False
                   				id_ = line.strip()
           				else:
                    				print('WARNING: An ID was found 
							without corresponding 
							sequence.', id_)
                    				id_ = line.strip()
            			elif id_ and not sequence:
               				sequence = line.strip()
            			elif id_ and sequence:
                   			score = line.strip()
                    			data[id_] = {'sequence':sequence,
                                	  		'score':score}
                    			sequence = False 
                 	   		id_ = False
    		return data

The data returns a dictionary with identification (id), sequence and scores for individual bases in the fastq file . Following this script, plotting was conducted as shown in figure 1. Here the distribution of the error rate is visualized for every individual base in the sequence, which illustrates the occurence of high and low score regions. In these high score regions bases have more chance to be falsly assigned and are thus less reliable for mutation identification. Likewise, mutations found in low score regions are more likely to be rightly assigned and thus true mutations. 

.. figure::  C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Fastq_files_qualityscore.png
   :scale:   70%
   :align:   center

   Figure 2: Distribution of error rate (quality score) for every individual base. 
   The mean of all the quality scores is 0.21.

While figure 1 gives a good overview on the importance of analyzing quality scores, one fastq file has little quantifiable value for identifying high score regions. To identify high score regions, 1139 fastq files from nanopore with data on chromosome 9 were likewise investigated for quality score distribution. Firstly, meanscores of the quality scores were plotted in figure 2. The meanscores of the files is a good first indication on the sequence quality and can be used to filter and select files for variant calling, where a low quality score is important. 

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Mean_distribution_of_all_quality_score.png
   :scale:  70%
   :align:  center

   Figure 3: Mean distribution of quality scores from all the fastq files.

Next, regions of bases were selected instead of single bases to be able to identify high quality score regions. Size and overlap of the chunks of sequences could be selected by the
following script::

	def split_overlap(iterable,size,overlap):
    		if size < 1 or overlap < 0:
        		raise ValueError('"size" must be an integer with
					 >= 1 while "overlap" must be >= 0')
    		result = []
    		while True:
        		if len(iterable) <= size:
            			result.append(iterable)
            			return result
        		else:
            			result.append(iterable[:size])
            			iterable = iterable[size-overlap:] 

Sequences were chunked to pieces of four to six in order to analyze the impact of different sized regions on base calling quality. Following, the mean of all the qualityscores of the same chunks of sequences were either plotted directly (figure 3 A - C) or indirectly after being devided in categories of high, medium and low qualityscore (figure 3 D - F). Categorizing was done after calculating mean of sequence, subsequently, categories were counted for each sequence. Categorizing was conducted to manipulate and increase data analysis. Parameters for categorizing were randomly selected and differentiate for each size, because with larger regions, the mean of the qualityscore get's more normalized and shift further towards medium, which have been accounted for by lowering high requirements and highering low requirements as following::

	def high_medium_low_scores(listed_scores, size):
    		group_score = []
    		for s in listed_scores:
        		if s >= (0.40-0.02*size):
           			group_score.append('High')
        		elif s <= (0.15+0.01*size):
            			group_score.append('Low')
        		else:
            			group_score.append('Medium')
    		return group_score

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\Source\\_static\\Fastq_gridplots.png
   :scale:  30%
   :align:  center

   Figure 4: **Quality score analysis with 6 senario's.** A - C) Meanscore for all combination in size (A = 4, B = 5, C = 6) for 1139 fastq files derived from nanopore sequencing of
   chromosome 9. D - F) Scores for regions have been categorized into high, medium and low for regions of same size as A to C. Next, the amount of times a region was called under a 
   certain category was counted and collected for the same data set. In these figures scores are set in percentage of total amount of times a region occurs in the data set.
   (Interactive figure at GridPlot_)

In table 1, highest and lowest five scoring sequence are highlighted. In conclusion, the highest scoring sequence has the biggest chance to have wrongly assigned bases in it's sequence.
In contrast, bases in lower scoring sequences are more likely to been good assigned and are therefor indeed the right base. These findings should be taken into account when investigating 
mutations, as a mutation found in for instance TTCC are more likely to be a real mutations than a mutation found in GCTT.

+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|  sequence | Score |  sequence | Score |  sequence | Score |  sequence | Score |  sequence | Score |  sequence | Score |
|     A     |       |     D     |   %   |     B     |       |     E     |   %   |     C     |       |     F     |   %   |
+===========+=======+===========+=======+===========+=======+===========+=======+===========+=======+===========+=======+
|   GCTT    | 0.364 |    GCTT   | 59.41 |   AGCTT   | 0.422 |   CCTTG   | 66.00 |   AGCTTT  | 0.501 |   TCATAC  | 91.52 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   CTTG    | 0.353 |    CTTG   | 58.04 |   GCCTT   | 0.405 |   CTTGC   | 65.52 |   TTCGCA  | 0.499 |   AGCCTT  | 90.00 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   TAAT    | 0.313 |    TAAT   | 46.72 |   GCTTG   | 0.393 |   CTTTA   | 65.00 |   GGGACG  | 0.489 |   CTTTAC  | 88.88 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   GTAG    | 0.298 |    GTAG   | 43.12 |   GCTTA   | 0.372 |   GTAGC   | 64.38 |   CCATGT  | 0.482 |   TAGCCA  | 87.50 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   TAGC    | 0.293 |    TAGC   | 42.61 |   ATTGA   | 0.367 |   CGGAG   | 63.16 |   GAATCT  | 0.466 |   TGCTAC  | 83.33 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|    ...    |       |    ...    |       |    ...    |       |    ...    |       |    ...    |       |    ...    |       |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   GGAT    | 0.136 |    GGTT   |  3.82 |   TTAAA   | 0.112 |   CGGGA   |  3.92 |   CCTAAT  | 0.058 |   TCCACT  |  1.33 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   CCCT    | 0.135 |    CCTC   |  3.64 |   GTCTT   | 0.104 |   CTCCT   |  3.88 |   TTCACA  | 0.054 |   TTATCC  |  1.23 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   GTTC    | 0.131 |    ATCC   |  3.53 |   TTGGA   | 0.100 |   CTCCA   |  2.93 |   TTTTTC  | 0.053 |   CCTCCT  |  1.18 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   CCTC    | 0.129 |    GATC   |  3.35 |   GGACC   | 0.098 |   CGATC   |  2.89 |   CCAATC  | 0.050 |   TCGGAT  |  1.05 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+
|   TTCC    | 0.128 |    CTCC   |  2.79 |   TTTTT   | 0.085 |   TCGGA   |  1.62 |   GGACGT  | 0.049 |   GGGACC  |  0.96 |
+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+-----------+-------+

   Table 1: Highest and lowest five scoring sequences. A - C) score is meanscore for all combinations in same size and data set as figure 3. 
   D - F) score is percentage of sequence in category high for all combinations in same size and data set as figure 3. 

**#Clustering**

Variant Call Format files
+++++++++++++++++++++++++
**Mutation distribution of single nucleotide polymorphisms**. Variant Call Format (VCF) files are text files containing data of single positions in the genome. In these files, variants
are formatted with the reference included. For sequenced sites, amount of reads found with mutation and reference are given. The dataset visualized here is derived from the cyclomics project, sequencing was preformed with nanopore and the data contains a sequence part from chromosome 17 (around 160 nucleotides) and a backbone, which is used for circulair pcr reaction. In total 1187 VCF files were used for variant calling. Here, VCF files are screened for single nucleotide polymorphism (SNP) occurence. Firsly, files were stripped of reported mutated bases, other data was discarded. As described earlier, every variant site has a number of reads that covers this site. These reads can be both coupled to the mutation and the reference. For example, on position 7577503 a SNP was found in 6 reads and 3 reads were coupled to the reference. While the amount of reads coupled to the mutation in contrast to the reads is important, here occurence of certain SNPs have been firsly investigated. In order to investigated the amount of SNPs in the files, VCF files were simallarly stripped as Fastq files and seperated by either sequence or backbone. Next, for the variants a parameter was set at a minimum of 25 percent of the reads that should be coupled to the mutant variant and visualized in figure 4:

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Combined_vcf_snp_analysis.png
   :scale:  70%
   :align:  center

   Figure 5: Distriution of SNPs in the sequence of 1187 VCF files. Parameter for variant identification was set at 25% of the reads to the variant. Variants are displayed as C -> T, meaning that T subsitutes C. A) Bar plot with single nucleotide polymorphisms occurence as percentage of whole. B) Heatmap from same variances with amount of occurences in the files

Both figures illustrate the common occurrence of G -> A mutation and to lesser extend due to C -> A. #what is more prompt to happen in normal biological system. #What could have happend here that made these mutations happen. #is there anyway we can change protocol to minimize these mutation frequenties. #difference between variances occuring in a single run (biased in that run) and variance occuring everytime by the machine (baised towards variance every run).

Next, SNPs were selected including 4 surrounding bases for heatmap analysis. Pandas was used to create a dataframe for the amount of times mutation occured to either A, T, C or G. This dataframe was then mapped with to a heatmap with reference sequence. #just as previous lenght can be changed for surrounding bases

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Variance_occurence_in_sequence_vcf.png
   :scale:  50%
   :align:  center

   Figure 5: Occurence of variance per reference sequence to different bases. In all the sequences the middle base is reported to be mutated in some of the vcf files. This mutation again has a parameter that is set at 25% of the reads atleast mutated. 

.. _GridPlot: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\gridplot.html

Script Tests
++++++++++++
Before scripts are ran over multiple files and directories, they should be checked for quality. In order to check a script for it's functionality, test scripts can be written. These testing scripts use the assert function, to identify if the set criteria are met.
As an example the earlier described parse_fasta_file_error is checked for it's quality with the following testing script::

	class TestDoneFastqParser:
    
    		def setup_method(self):
        		sequence_file = 'C:/Users/Douwe/Documents/Python/test_cases/test_fastq2.done_fastq'
        		self.data = dl.parse_fasta_file_error(sequence_file)
        		id_ = list(self.data.keys())[0]
        		self.score = self.data[id_]['score']

    		def check_valid_DNA_sequence(self, s):
        		for l in set(s.upper()):
            			if not l in 'ACTGN':
                			return False
        		return True
        
    		def test_has_id(self):
        		for id in '@':
            			assert id in list(self.data.keys())[0]
           
    		def test_sequence_correct(self):
        		for k, v in self.data.items():
            			assert self.check_valid_DNA_sequence(v['sequence']) == True
            
    		def test_score_correct(self):
        		for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            			assert letter not in self.score 

The class function is used to define which script is going to be checked for quality. Firstly the script is setup with a test file, this file is designed to identify flaws in the script. In other words, it consists off alot of errors which the script should not pickup. Next, multiple assertions are made, such as the assertion that letters in sequence can only consist of A, C, T, G and N. Also score should consist of characters and not involve any letters.
While this is an example of a test script, multiple scripts have been investigated for quality as described in the supplementairy.

References 
----------
(Ask how to link ref in text)

1: Gerlinger, M. et al. Intratumor Heterogeneity and Branched Evolution Revealed by Multiregion Sequencing. N. Engl. J. Med. 366, 883-892 (2012).
2: Choi, J.-J., Reich, C. F., Pisetsky, D. S. & Pisetsky, D. S. The role of macrophages in the in vitro generation of extracellular DNA from apoptotic and necrotic cells. Immunology 115, 55-62 (2005).
3: Beaver, J. A. et al. Detection of Cancer DNA in Plasma of Patients with Early-Stage Breast Cancer. Clin. Cancer Res. 20, 2643-2650 (2014).
4: Paper of Kloosterman group (not yet released)
5: Newman, A. M., Lovejoy, A. F., Klass, D. M., Kurtz, D. M., Chabon, J. J., Scherer, F., A. A. (2016). Integrated digital error suppression for improved detection of circulating tumor DNA. Nature Biotechnology, 34(5), 547-555.


.. toctree::
   :maxdepth: 2
   :caption: Contents: