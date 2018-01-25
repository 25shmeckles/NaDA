Results and discussion
----------------------
Fastq files
+++++++++++
**Quality score distribution**. Fastq are text-based files derived from sequencers such as Illumina and Nanopore that contain a biological sequence and for every base a quality score.
this quality score represents the sequencers ability to call a base as true, in other words, the score says something about how sure the sequencer is that the right base is called. 
if the score for the base is low, it's most likely the right base, while if the score is high the machine is unsure if it's the real base. Analyzing quality scores on sequences
is important to be able to understand which regions of bases is difficult to get base called. More importantly, a mutation called in a region with low quality score is therefore more likely to be a real mutation in contrast to a mutation that occurs through wrong base calling in area's with high quality scores. Clearifying which regions have high scores is important for achieving tests with high specificity, as these regions will not be selected for variant calling. Likewise, regions with low score are good targets for identifying mutations and thus are able to increase tests sensitivity. In order to analyse fastq files the following function was written to strip the files:: 

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

The data returns a dictionary with identification (id), sequence and scores for individual bases in the fastq file. Following this script, plotting was conducted as shown in figure 1. Here the distribution of the error rate is visualized for every individual base in the sequence, which illustrates the occurence of high and low score regions. In these high score regions bases have more chance to be falsly assigned and are thus less reliable for mutation identification. Likewise, mutations found in low score regions are more likely to be rightly assigned and thus true mutations. 

.. figure::  C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Fastq_files_qualityscore.png
   :scale:   70%
   :align:   center

   Figure 2: Distribution of error rate (quality score) for every individual base. Here every 10th base is visualized on the X-axis to clearify the data. The mean of all the quality scores is 0.21.

While figure 1 gives a good overview on the importance of analyzing quality scores, one fastq file has little quantifiable value for identifying high score regions. To identify high score regions, 1139 fastq files from nanopore with data on chromosome 9 were likewise investigated for quality score distribution. Firstly, meanscores of the quality scores were plotted in figure 2. The meanscores of the files is a good first indication on the sequence quality and can be used to filter and select files for variant calling, where a low mean quality score is important. 

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Mean_distribution_of_all_quality_score.png
   :scale:  100%
   :align:  center

   Figure 3: Mean distribution of quality scores from all the fastq files.

Next, regions of bases were selected instead of single bases to be able to identify high quality score regions. Size and overlap of the chunks of sequences could be selected by the
following function::

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

Sequences were chunked to pieces of four (tetrameer) to six (hexameer) in order to analyze the impact of different sized regions on base calling quality.
#why this size
Following, the mean of all the qualityscores of the same chunks of sequences were either plotted directly (figure 3 A - C) or indirectly after being devided in categories of high, medium and low qualityscore (figure 3 D - F). Categorizing was done after calculating mean of sequence, subsequently, categories were counted for each sequence. Categorizing was conducted to manipulate and increase data analysis. Parameters for categorizing were randomly selected and differentiate for each size, because with larger regions, the mean of the qualityscore get's more normalized and shift further towards medium, which have been accounted for by lowering high requirements and highering low requirements as following::

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

   Figure 4: **Quality score analysis with 6 senario's.** A - C) Meanscore for all combination in size (A = 4, B = 5, C = 6) for 1139 fastq files derived from nanopore sequencing of chromosome 9. D - F) Scores for regions have been categorized into high, medium and low for regions of same size as A to C. Next, the amount of times a region was called under a certain category was counted and collected for the same data set. In these figures scores are set in percentage of total amount of times a region occurs in the data set.(Interactive figure at GridPlot_)

.. _GridPlot: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\gridplot.html

In table 1, highest and lowest five scoring sequence are highlighted. In conclusion, the highest scoring sequence has the biggest chance to have wrongly assigned bases in it's sequence.
In contrast, bases in lower scoring sequences are more likely to been good assigned and are therefor indeed the right base. These findings should be taken into account when investigating 
mutations, as a mutation found in for instance TTCC is more likely to be a real mutations than a mutation found in GCTT.

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
|   ...     |       |    ...    |       |    ...    |       |    ...    |       |    ...    |       |    ...    |       |
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

**Clustering**. Another way of visualizing the qualityscore in fastq files is by using clustering. Clustering is a method in which data point get coupled in groups (clusters) by a certrain geometry. Here K-Means is used for clustering, which makes clusters based on the the following algorithm:

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\Source\\_static\\K-Means.png
   :scale:  70%
   :align:  center

#explain what this algorithm does

in figure 5 three clusters are formed in which the yellow one represents sequences with often reported high score and few times reported low scores. In this cluster, alterations are more likely to be falsly assigned. Furthermore, in blue cluster, alterations are more likely to be rightfully assigned. Clustering of data can provide for a more clearer view on which sequences to include and exclude for mutation calling.

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\Source\\_static\\clusterplot.png
   :scale:  50%
   :align:  center

   figure 5: **Clustering of data from 1139 fastq files.** Sequences of 5 bases are measured for qualityscore and reported in high or low score. Here the percentage of times sequence is reported in high and low score is visualized. Following, clustering was conducted using K-means, seperating three clusters. 

Together fastq data suggest that qualityscore is important in identifying regions which are promosing for mutation calling and which regions should be avoided. As described earlier, regions with a high qualityscore should be avoided while looking for mutation. In contrast, low qualityscore region have potential for identifying mutations in cfDNA.

Importantly, the qualityscores of sequences can differ on the method being used. In this case our method involves rolling circle amplifcation and nanopore sequencing of cfDNA. In order to make a sensitive data filter, a big database of healthy cfDNA should investigated on qualityscore for sequences. Therefore, the filter can exclude and include regions with respectivilly high and low qualityscores. Furthermore, qualityscores of sequences can differ on every run, causing some sequences to have higher or lower scores. For this discrepancy should also be accounted in the data filter. A possiblity is to either include healthy cfDNA into every run or compare backbone sequence to identify run specific sequence qualityscore differences.

Here only run specific sequence qualityscores have been investigated. In order to visualize high and low score sequences specific for our method, data analysis should be conducted on multiple runs. Analysis of multiple runs can be simultaneously conducted using the High-Performance Computing (HPC) facility in the UMC, which will be done in the following segment.
On the HPC multiple ctDNA datasets derived from cyclomics are available for analysis, here the focus goes towards the rolling circle amplification p53 mutated and wild-type(WT) datasets. Firstly, a fastq_script_ was written to achieve similair data analysis and visualization as described above. Minor visualization updates were conducted to improve data visibility. Both datasets are seperated in equaly sized chunks, around 4000 files each, and analysed as individual chunks to increase script parallelization, thus increasing speed. For all files meanscores were calculated and visualized in figure 6. This boxplot clearly visualizes the lack of consistency between qualityscores in the same sequence run. Therefor, this could give an indication that qualityscores have limited value for developing a data filter. 

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\Source\\_static\\RCAxMUT_WT_boxplot.png
   :scale:  70%
   :align:  center

   Figure 6: **Boxplot of meanscore from several chunks in p53 mutated and WT dataset.** 

.. _fastq_script: https://github.com/DouweSpaanderman/NaDA/blob/master/Scripts/fastq_qualityscore_analyser.py

While meanscores give an indication on qualityscore analysis, both qualityscore plotting and clustering is yet to determine if high and low score region excist and percist in multiple chunks and datasets. For every chunk derived from a dataset, sequences have been analysed and visualized as tetrameer, pentameer and hexameer. Here, tetrameers of wild-type chunk 0 to 3 have been plotted as shown in figure 7.

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\Source\\_static\\Fastq_gridplot_WT.png
   :scale:  70%
   :align:  center

   Figure 7: **tetrameer sequence analysis for chunks zero to three of the p53 wild-type database.** Figures illustrates the mean qualityscore for each tetrameer possible in one chunk. A) chunk0. B) chunk1. C) chunk2. D) chunk3. Interactive figure can be found here and also visualizes data analysis when devided into high, medium and low group.(WT_chunk0_, WT_chunk1_, WT_chunk2_ and WT_chunk3_)

.. _WT_chunk0: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk0_4.0_3.0_score_plotting.html
.. _WT_chunk1: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk1_4.0_3.0_score_plotting.html
.. _WT_chunk2: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk2_4.0_3.0_score_plotting.html
.. _WT_chunk3: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk3_4.0_3.0_score_plotting.html

Similarly to the boxplot, their seems to be a lack of consistency between chunks as high reported tetrameers differ heavely between these chunks. Similair results are visuable for bigger sized chunks(supplementary_1_) and chunks derived from p53 mutant dataset(supplementary_2_). These datasets show that their is yet to be proven for a correlation between qualityscores and specific regions or chunks. However, clustering could clarify for the occurence of high qualityscore regions by better identification of these regions. In order to cluster data derived from dataset chunks, the same algorithm is used as described above. In figure 8 clustering of chunks 0 to 3 has been visualized.

.. _supplementary_1: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\build\\html\\Supplementary.html
.. _supplementary_2: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\build\\html\\Supplementary.html

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\Source\\_static\\Fastq_gridplot_WT_cluster.png
   :scale:  70%
   :align:  center

   Figure 8: **Clustering of hexameer sequence for chunks zero to three of the p53 wild-type database.** A) chunk0. B) chunk1. C) chunk2. D) chunk3. Interactive cluster plot can be found here. (WT_chunk0_cluster_, WT_chunk1_cluster_, WT_chunk2_cluster_ and WT_chunk3_cluster_)

.. _WT_chunk0_cluster: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk0_6.0_5.0_score_clustering.html
.. _WT_chunk1_cluster: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk1_6.0_5.0_score_clustering.html
.. _WT_chunk2_cluster: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk2_6.0_5.0_score_clustering.html
.. _WT_chunk3_cluster: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxWT_chunk3_6.0_5.0_score_clustering.html

#clustering figures

#conclusion and what's next with fastq

Variant Call Format files
+++++++++++++++++++++++++
**Mutation distribution of single nucleotide polymorphisms**. Variant Call Format (VCF) files are text files containing data of single positions in the genome. In these files, variants
are formatted with the reference included. For sequenced sites, amount of reads found with mutation and reference are given. The dataset visualized here is derived from the cyclomics project, sequencing was preformed with nanopore and the data contains a part of the p53 gen on chromosome 17 (around 160 nucleotides) and a backbone, which is used for circulair pcr reaction. In total 1187 VCF files were used for variant calling. Here, VCF files are screened for single nucleotide polymorphism (SNP) occurence. Firstly, files were stripped of reported mutated bases, other data was discarded. As described earlier, every variant site has a number of reads that covers this site. These reads can be both coupled to the mutation and the reference. For example, on position 7577503 a SNP was found in 6 reads and 3 reads were coupled to the reference. While the amount of reads coupled to the mutation in contrast to the reads is important, here occurence of certain SNPs have been firsly investigated. In order to investigated the amount of SNPs in the files, VCF files were simallarly stripped as Fastq files and seperated by either sequence or backbone. Next, for the variants a parameter was set at a minimum of 25 percent of the reads that should be coupled to the mutant variant and visualized in figure 4:

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Combined_vcf_snp_analysis.png
   :scale:  70%
   :align:  center

   Figure 6: Distriution of SNPs in the sequence of the p53 gen for 1187 VCF files. Parameter for variant identification was set at 25% of the reads to the variant. Variants are displayed as C > T, meaning that T subsitutes C. A) Bar plot with single nucleotide polymorphisms occurence as percentage of whole. B) Heatmap from same variances with amount of occurences in the files

Both figures illustrate the common occurrence of G > A mutation and to lesser extend due to C > A. The prevalance of these SNPs in contrast to other alterations are a strong indication that these alterations are caused by a non-biological mechanism, which can be errors in the rolling circle amplification, library preparation and sequencing of the ctDNA. In literature, cytosine deamination has been described to increase C:G > T:A noise levels (6). Also, less occurring alteration C > A has been reported to be caused by oxidative DNA damage during sample preparation(7). Both these types of alterations can be a result of polymerase-induced errors. Possible suggested methods to suppress these errors are adding DNA repair mechanisms upon polymerase chain reaction (PCR) and lowering heat. However, an in silico approach to polish background noise can also be devised. 

Next, SNPs were selected in a triplet for heatmap analysis. Pandas was used to create a dataframe for the amount of times mutation occured to either A, T, C or G. This dataframe was then mapped to a heatmap with reference sequence. Just as in previous figures, lenght of the surrounding bases can be changed to give a wider variety of information. This gave more information about base combinations with high alteration affinity, such as ACGCA to ACACA. 

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Variance_occurence_in_sequence_vcf_3.png
   :scale:  70%
   :align:  center

   Figure 7: Occurence of variance per reference sequence to different bases. In all the sequences the middle base is reported to be mutated in some of the vcf files. This mutation again has a parameter that is set at 25% of the reads atleast mutated. 

Identifying high variance regions in both healthy cfDNA and ctDNA is important for constructing a data filter. #something

Furthermore, just as with the fastq files, variances can be seperated between alterations specific for a run and alterations specific for the method being used. For instance, CTC > A could be an alteration that is specifically highly mutated in a particularly run, while CGC > A occurs often in every run with this method of rolling circle amplification and nanopore sequencing. Therefore, filtering should be able to account for both run specific and method specific alterations. In the same manner, high database of healthy cfDNA could accomplish a method specific filter and adding healthy cfDNA into every run a specific alterations filter. Also more covenient, backbone data could be used to identify run specific errors as the backbone doesn't change between runs and should thus never contain alterations.

Further analysis on VCF files was conducted on the HPC system. In order to conduct VCF analysis on the HPC, a vcf analyis script_ was written. This script analysis mutation occurence as a percentage of amount of times sequence has been reported in the dataset. This script was run over multiple Datasets available from cyclomics, which were p53 wild-type(WT), mutated(MUT), 1% mutated and 10% mutated. Vcf script were also directly visualized the dataset in heatmap and barplot. Here, the p53 mutated dataset has been visualized in figure x. As expected, mutation occurence (GGATA > T) can be clearly visualised in this dataset. Furthermore, WT analysis showed a compleet absence of this mutation (supplementary_a_) and the other two dataset are in concordance with the percentage mutated (supplementary_b_). 

Extraordinarily, apart from the expected mutation, in all four datasets other pentameers seem to be mutated aswell. Especially, CAACC is reported to be highly mutated (around 30%) for all the datasets. This could indicate either the occurence of other mutations in the dataset or the identification of pentameers which cause error's throughout our workflow.

.. figure:: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\Variance_occurence_in_sequence_vcf_3.png
   :scale:  70%
   :align:  center

   Figure x: Occurence of alterations as a percentages of the amount of times sequence has been reported in the dataset. In all the sequences the middle base has been reported to be mutated in some of the vcf files. Dataset used here is p53 mutated. Interactive figure can be found here(MUT_heatmap_)

.. _supplementary_a: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\build\\html\\Supplementary.html
.. _supplementary_b: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\build\\html\\Supplementary.html
.. _script: https://github.com/DouweSpaanderman/NaDA/blob/master/Scripts/vcf_snp_variant_analyser.py
.. _MUT_heatmap: C:\\Users\\Douwe\\Documents\\GitHub\\NaDA\\Documentation\\source\\_static\\RCAxMUT_insert_5.0_heatmap_sequences.html

Importantly, the identified alterations could persist anywhere the pentameer is located in the sequence. Furthermore, if a alteration occurs in one specific position it is more prone to be an actual mutation. In contrast, alterations occuring in multiple location with the same pentameer could indicate to a systematic problem with pcr or sequencing causing an falsly identified mutation. In order to visualize the location of the mutation, both a heatmap covering single nucleotides instead of pentameers was constructed (figure x) and another script was prepaired to compare pentameer locations in the sequence(link of excel). 



Script Tests
++++++++++++
Before scripts are run over multiple files and directories, they should be checked for quality. In order to check a script for it's functionality, test scripts can be written. These testing scripts use the assert function to identify if the set criteria are met.
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

The class function is used to define which script is going to be checked for quality. Firstly the script is setup with a test file, this file is designed to identify flaws in the script. In other words, it consists off alot of errors which the script should not pickup. Next, multiple assertions are made, such as the assertion that letters in sequence can only consist of A, C, T, G and N. Also score should consist of characters and not involve any letters. While this is an example of a test script, multiple scripts have been investigated for quality as described in the supplementairy.