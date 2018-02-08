Results and discussion
----------------------
Fastq files
+++++++++++
**Quality score distribution**. Fastq are text-based files derived from sequencers such as Illumina and Nanopore that contain a biological sequence and for every base a quality score.
this quality score represents the sequencers ability to call a base as true, in other words, the score says something about how sure the sequencer is that the right base is called. 
if the score for the base is low, it's most likely the right base, while if the score is high the machine is unsure if it's the real base. Analyzing quality scores on sequences
is important to be able to understand which regions of bases are difficult to get base called. More importantly, a mutation called in a region with low quality score is therefore more likely to be a real mutation in contrast to a mutation that occurs through wrong base calling in area's with high quality scores. Clarifying which regions have high scores is important for achieving tests with high specificity, as these regions will not be selected for variant calling. Likewise, regions with low score are good targets for identifying mutations and thus are able to increase test sensitivity. In order to analyse fastq files the following function was written to strip the files:: 

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

The data returns a dictionary with identification (id), sequence and scores for individual bases in the fastq file. Following this script, the plotting was conducted as shown in figure 1. Here the distribution of the error rate is visualized for every individual base in the sequence, which illustrates the occurrence of high and low score regions. In these high score regions, bases have more chance to be falsely assigned and are thus less reliable for mutation identification. Likewise, mutations found in low score regions are more likely to be rightly assigned and thus true mutations. 

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Fastq_files_qualityscore.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 2: **Distribution of error rate (quality score) for every individual base.** Here every 10th base is visualized on the X-axis to clarify the data. The mean of all the quality scores is 0.21.

While figure 1 gives a good overview on the importance of analyzing quality scores, one fastq file has little quantifiable value for identifying high score regions. To identify high score regions, 1139 fastq files from Nanopore with reads on chromosome 9 were likewise investigated for quality score distribution. Firstly, mean scores of the quality scores were plotted in figure 2. The mean scores of the files are a good first indication of the sequence quality and can be used to filter and select files for variant calling, where a low mean quality score is important. 

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Mean_distribution_of_all_quality_score.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 3: **Mean distribution of quality scores**. Mean scores from 1139 fastq files were calculated and visualised in a violin plot.

**Region quality score distribution**. Next, regions of bases were selected instead of single bases to be able to identify high quality score regions. Size and overlap of the chunks of sequences could be selected by the following function::

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

Sequences were chunked to pieces of four (tetrameer) to six (hexameer) in order to analyze the impact of different sized regions on base calling quality. Following, the mean of all the quality scores of the same chunks of sequences were either plotted directly (figure 3 A - C) or indirectly after being divided in categories of high, medium and low quality score (figure 3 D - F). Categorizing was done before calculating mean of a sequence, subsequently, categories were counted for each sequence. Categorizing was conducted to manipulate and increase data analysis. Parameters for categorizing were randomly selected and differentiate for each size, because with larger regions, the mean of the quality score get's more normalized and shift further towards medium, which have been accounted for by lowering high requirements and increasing low requirements as following::

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

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Fastq_gridplots.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 4: **Quality score analysis in 6 scenarios.** A - C) Mean score for all combinations in chunk sizes (A = 4, B = 5, C = 6) for 1139 fastq files derived from nanopore sequencing of chromosome 9. D - F) Scores for regions have been categorized into high, medium and low for regions of the same size as A to C. Next, the amount of times a region was called under a certain category was counted and collected for the same data set. In these figures, scores are set in percentage of the total amount of times a region occurs in the data set.(Interactive figure at GridPlot_)

.. _GridPlot: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/gridplot.html

In table 1, highest and lowest five scoring sequences are highlighted. In conclusion, the highest scoring sequence has the biggest chance to have wrongly assigned bases in its sequence. In contrast, bases in lower scoring sequences have been most likely rightfully assigned and are therefore indeed the right base. These findings should be taken into account when investigating mutations, as a mutation found in, for instance TTCC, is more likely to be a real mutation than a alterations found in GCTT.

+-----------+-------+-----------+-------+-----------+-------+
|     C     | Score |     B     | Score |     C     | Score |
|           |       |           |       |           |       |
+===========+=======+===========+=======+===========+=======+
|   GCTT    | 0.364 |   AGCTT   | 0.422 |   AGCTTT  | 0.501 |
+-----------+-------+-----------+-------+-----------+-------+
|   CTTG    | 0.353 |   GCCTT   | 0.405 |   TTCGCA  | 0.499 |
+-----------+-------+-----------+-------+-----------+-------+
|   TAAT    | 0.313 |   GCTTG   | 0.393 |   GGGACG  | 0.489 |
+-----------+-------+-----------+-------+-----------+-------+
|   GTAG    | 0.298 |   GCTTA   | 0.372 |   CCATGT  | 0.482 |
+-----------+-------+-----------+-------+-----------+-------+
|   TAGC    | 0.293 |   ATTGA   | 0.367 |   GAATCT  | 0.466 |
+-----------+-------+-----------+-------+-----------+-------+
|   ...     |       |    ...    |       |    ...    |       |
+-----------+-------+-----------+-------+-----------+-------+
|   GGAT    | 0.136 |   TTAAA   | 0.112 |   CCTAAT  | 0.058 |
+-----------+-------+-----------+-------+-----------+-------+
|   CCCT    | 0.135 |   GTCTT   | 0.104 |   TTCACA  | 0.054 |
+-----------+-------+-----------+-------+-----------+-------+
|   GTTC    | 0.131 |   TTGGA   | 0.100 |   TTTTTC  | 0.053 |
+-----------+-------+-----------+-------+-----------+-------+
|   CCTC    | 0.129 |   GGACC   | 0.098 |   CCAATC  | 0.050 |
+-----------+-------+-----------+-------+-----------+-------+
|   TTCC    | 0.128 |   TTTTT   | 0.085 |   GGACGT  | 0.049 |
+-----------+-------+-----------+-------+-----------+-------+

   Table 1: **Highest and lowest five scoring sequences**. A - C) score is mean score for all combinations of same size and data set as figure 3.

**Clustering**. Another way of visualizing the quality score in fastq files is by using clustering. Clustering is a method in which data points get coupled in groups (clusters) by a certain geometry. Here K-Means is used for clustering, which makes clusters based on the distance between data points and a selected centroid. The centroid is the mean of a cluster and is defined by a trial and error process. This process is repeated until centroids are selected, which happens when the within-cluster sum of squares is minimised.

In figure 5 three clusters are formed in which the yellow one represents sequences with often reported high score and few times reported low scores. In this cluster, alterations are more likely to be falsely assigned. Furthermore, in blue cluster, alterations are more likely to be rightfully assigned. Clustering of data can provide for a more clearer view on which sequences to include and exclude for mutation calling.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/clusterplot.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 5: **Clustering of data from 1139 fastq files.** Sequences of 5 bases (pentameer) are measured for quality score and reported in high or low score. Here the percentage of amount a sequence is reported in high and low score is visualized. Following, clustering was conducted using K-means, separating three clusters. 

Together fastq data suggest that quality score is important in identifying regions which are promising for mutation calling and which regions should be avoided. As described earlier, regions with a high quality score should be avoided while looking for mutation. In contrast, low quality score regions have potential for identifying mutations in cfDNA.

Importantly, the quality scores of sequences can differ on the method being used. In this case our method involves rolling circle amplification and nanopore sequencing of cfDNA. In order to make a sensitive data filter, a big database of healthy cfDNA should be investigated on quality score for sequences. Therefore, the filter can exclude and include regions with respectively high and low quality scores. Furthermore, quality scores of sequences can differ on every run, causing some sequences to have higher or lower scores. For this discrepancy should also be accounted in the data filter. A possibility is to either include healthy cfDNA into every run or compare backbone sequences to identify run specific sequence quality score differences.

**p53 wild-type and mutant dataset analysis**. So far, only run specific sequence quality scores have been investigated. In order to visualize high and low score sequences specific for our method, data analysis should be conducted on multiple runs. Analysis of multiple runs can be simultaneously conducted using the High-Performance Computing (HPC) facility in the UMC, which will be done in the following segment.

On the HPC multiple ctDNA datasets derived from cyclomics are available for analysis, here the focus goes towards the rolling circle amplification p53 mutated and wild-type(WT) datasets. Firstly, a fastq_script_ was written to achieve similar data analysis and visualization as described above. Minor visualization updates were conducted to improve data visibility. Both datasets are separated in equally sized chunks, around 4000 files each, and analysed as individual chunks to increase script parallelization, thus increasing speed. For all files mean scores were calculated and visualized in figure 6.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/RCAxMUT_WT_boxplot.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 6: **Boxplot of mean score from p53 mutated and WT dataset.** For both datasets chunk 0 to 9 have been visualized. Chunk 10 to 24 were excluded, but contained similar results.

.. _fastq_script: https://github.com/DouweSpaanderman/NaDA/blob/master/Scripts/fastq_qualityscore_analyser.py

This boxplot clearly visualizes the lack of consistency between quality scores in the same sequence run. Therefor, this could give an indication that quality scores have limited value for developing a data filter. 

While mean scores give an indication on quality score analysis, both quality score plotting and clustering is yet to determine if high and low score region exist and persist in multiple chunks and datasets. For every chunk derived from a dataset, sequences have been analysed and visualized as tetrameer, pentameer and hexameer. Here, tetrameers of wild-type chunk 0 to 3 have been plotted as shown in figure 7.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Fastq_gridplot_WT.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 7: **tetrameer sequence analysis for chunks zero to three of the p53 wild-type database.** Figures illustrates the mean qualityscore for each tetrameer possible in one chunk. A) chunk 0. B) chunk 1. C) chunk 2. D) chunk 3. Interactive figure can be found here and also visualizes data analysis when divided into high, medium and low group.(WT_chunk0_, WT_chunk1_, WT_chunk2_ and WT_chunk3_)

.. _WT_chunk0: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk0_4.0_3.0_score_plotting.html
.. _WT_chunk1: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk1_4.0_3.0_score_plotting.html
.. _WT_chunk2: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk2_4.0_3.0_score_plotting.html
.. _WT_chunk3: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk3_4.0_3.0_score_plotting.html

Similarly to the boxplot, there seems to be a lack of consistency between chunks as high reported tetrameers differ heavily between these chunks. Similar results are visible for bigger sized chunks(supplementary_1_) and chunks derived from the p53 mutant dataset(supplementary_2_). These datasets show that their is yet to be proven that a correlation between quality scores and specific regions or chunks exists. However, clustering could clarify the occurrence of high quality score regions by better identification of these regions. In order to cluster data derived from dataset chunks, the same algorithm is used as described above. In figure 8 clustering of chunks 0 to 3 from p53 WT has been visualized.

.. _supplementary_1: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html
.. _supplementary_2: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Fastq_gridplot_WT_cluster.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 8: **Clustering of hexameer sequence for chunks zero to three of the p53 wild-type database.** Data points are visualised as percentage reported in high (y-axis) and low score(x-axis). A) chunk 0. B) chunk 1. C) chunk 2. D) chunk 3. Interactive cluster plot can be found here. (WT_chunk0_cluster_, WT_chunk1_cluster_, WT_chunk2_cluster_ and WT_chunk3_cluster_)

.. _WT_chunk0_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk0_6.0_5.0_score_clustering.html
.. _WT_chunk1_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk1_6.0_5.0_score_clustering.html
.. _WT_chunk2_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk2_6.0_5.0_score_clustering.html
.. _WT_chunk3_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk3_6.0_5.0_score_clustering.html

Clustering is able to identify regions that have both been reported often as high score and few times as low scores. However, between chunks there is a huge discrepancy in quality scores. Chunk 0 and 2 have an overall much lower quality score in comparison with chunk 1 and 2, which was also identified with mean quality scores in figure 6. This big difference in overall quality of the chunks seems to indicate that quality score can't be used for constructing a data filter. Nevertheless, if regions are always present in the same cluster between chunks, high and low quality score regions could still be identified. Additionally, quality score mean shouldn't have to influence score from a single region in comparison to other regions in the same chunk.

Identification of these regions could help data filtering for mutation, as bases in regions with high quality score are less likely to be rightfully assigned and the other way around for low quality score region. Therefore respectively, these regions could be excluded or included in variant calling. Currently, quality score analysis shows a lot of inconsistency between chunks of the same dataset. Thus, it is yet impossible to conclude any regions that have either a high or low quality score. Therefor, quality score has currently no application in creating a data filter. with all this in mind, while quality score shows limited possibilities, regions should still be compared between chunks for in which cluster they are reported. This could clarify if there is a correlation between regions and quality scores.

Variant Call Format files
+++++++++++++++++++++++++
**Mutation distribution of single nucleotide polymorphisms**. Variant Call Format (VCF) files are text files containing data of single positions in the genome. In these files, variants are formatted with the reference included. For sequenced sites, amount of reads found with the mutation and reference are given. The dataset visualized here is derived from the Cyclomics project, sequencing was preformed with Nanopore and the data contains a part of the p53 gen on chromosome 17, around 160 nucleotides and a backbone, which is used for circular pcr reaction. In total 1187 VCF files were used for variant calling. Here, VCF files are screened for single nucleotide polymorphism (SNP) occurrence. Firstly, files were stripped of reported mutated bases, other data was discarded. As described earlier, every variant site has a number of reads that covers this site. These reads can be both coupled to the mutation and the reference. For example, on position 7577503 a SNP was found in 6 reads and 3 reads were coupled to the reference. While the amount of reads coupled to the mutation in contrast to the total amount of reads is important, here the occurrence of certain SNPs have been firstly investigated. Analysing SNPs in ctDNA could help identify real mutation present in the tumour. Likewise, analysing cfDNA for SNPs can be useful for lowering background noise, by identifying passenger mutation, asymmetric DNA errors, PCR errors and sequencing errors. 

In order to investigate the amount of SNPs in the files, VCF files were similarly stripped as Fastq files and separated by either sequence or backbone. Next, for the alterations a parameter was set at a minimum of 25 percent of the reads that should be coupled to the mutant variant and visualized in figure 9:

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Combined_vcf_snp_analysis.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 9: **Distriution of SNPs in the sequence of the p53 gen for 1187 VCF files.** Parameter for variant identification was set at 25% of the reads to the variant. Variants are displayed as C > T, meaning that T substitutes C. A) Bar plot with single nucleotide polymorphisms occurrence as percentage of whole. B) Heatmap from same variances with amount of occurrences in the files

Both figures illustrate the common occurrence of G > A mutation and to lesser extend C > A. The prevalence of these SNPs in contrast to other alterations are a strong indication that these alterations are caused by a non-biological mechanism, which can be errors in the rolling circle amplification, library preparation and sequencing of the ctDNA. In literature, cytosine deamination has been described to increase C:G > T:A noise levels (6). Also, less occurring alteration C > A has been reported to be caused by oxidative DNA damage during sample preparation(7). Both these types of alterations can be a result of polymerase-induced errors. Possible suggested methods to suppress these errors are adding DNA repair mechanisms upon polymerase chain reaction (PCR) and lowering heat. However, if results are similar in other databases, an in silico approach to polish background noise can also be devised. 

**Region SNP analysis**. Next, trimeers and pentameers were selected in which the middle base was reported to contain a SNPs in some of the files for heatmap analysis. Creating chunks for SNP analysis instead of single base analysis was conducted to visualize sequences that were more likely to contain SNPs. Pandas library was used to create a dataframe for the amount of times mutation occured to either A, T, C or G. This dataframe was then mapped to a heatmap with reference sequence. Just as in previous figures, lenght of the surrounding bases can be changed to give a wider variety of information. This gave more information about base combinations with high alteration affinity, such as ACGCA to ACACA. 

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Variance_occurence_in_sequence_vcf_3.PNG
   :width: 700px
   :height: 450px
   :align: center

   Figure 10: **Occurence of variance per reference sequence to different bases.** In all the sequences the middle base is reported to be mutated in some of the vcf files. This mutation again has a parameter that is set at 25% of the reads atleast mutated. 

Identifying high variance chunks in healthy cfDNA is important for understanding where errors arise in the sequence and thus could help polishing background noise by deselecting these chunks for variant calling analysis.

Furthermore, just as with the fastq files, variances can be separated between alterations specific for a run and alterations specific for the method being used. For instance, CTC > A could be an alteration that is specifically highly mutated in a particular run, while CGC > A occurs often in every run with this method of rolling circle amplification and Nanopore sequencing. Therefore, filtering should be able to account for both run specific and method specific alterations. In the same manner, a large database of healthy cfDNA could accomplish a method specific filter and adding healthy cfDNA into every run a specific alteration filter. Also more convenient, backbone data could be used to identify run specific errors as the backbone doesn't change between runs and should thus never contain alterations.

**p53 wild-type and mutant dataset analysis**. Further analysis on VCF files was conducted on the HPC system. In order to conduct VCF analysis on the HPC, a vcf analysis script_ was written. This script analysis mutation occurrence in insert and backbone similar as described above. Subsequently, visualization was updated, presenting sequences as a percentage of times it has been reported to contain a SNPs to how many times it was reported in total for the whole dataset. Furthermore, SNPs were also analysed as a single location rather than a chunk of three to five bases.

This script was run over multiple datasets available from Cyclomics, which were p53 wild-type(WT), mutated(MUT), 1% mutated and 10% mutated. Firstly, SNPs were investigated for their individual positions. As visualized in figure 11, using this single position heatmap, a clear mutation, position 7578265 A > T, can be found in the p53 mutated dataset. 

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/single_base_insert_MUT.PNG
   :width: 700px
   :height: 450px
   :align: center

   Figure 11: **Single position SNPs analysis for p53 mutated database**. Occurrence of SNP is visualized as a percentage of times it has been reported with a SNP to the total amount of times it has been reported in the dataset. Interactive figure can be found here(single_base_insert_MUT_)

.. _single_base_insert_MUT: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxMUT_single_base_insert_1_heatmap_sequences.html

Similar to figure 11, wild-type, 1% and 10% have been analysed for position specific SNPs (supplementary_3_). As expected, wild-type shows no occurrence of the specific mutation and the other datasets are in concordance with the percentages of reads that should contain the SNP. Therefor, this script_ is able to identify true mutations in this p53 dataset.

.. _supplementary_3: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html

However, while the expected alteration is found in our datasets, for those datasets other alteration were also identified. These alterations were not in line with our expectations and could have a number of explanations. Firstly, it could be a real passenger mutation which is not normally present in the reference. Secondly, these alterations could be caused by protocol errors, such as pcr errors and sequencing errors. For both these alterations could be accounted for by excluding certain positions in the data filter. Nevertheless, this could only help to construct a data filter specific for p53 databases as specific positions get excluded and our aim is to construct a data filter that can be applied to multiple gene databases. 

In order to construct a data filter usable for multiple genes, regions rather than single positions were investigated for SNPs. Region selection could help identifying SNPs due to protocol errors and help exclude them from analysis. First, chunks of trimeers and pentameers were similarly selected as described in region SNP analysis. Subsequently, these chunks were plotted in a heatmap visualized in figure 12. Mutation occurrence (GGATA > T) could again be clearly visualised in this dataset. Furthermore, WT analysis showed a complete absence of this mutation (supplementary_4_) and the other two datasets are in concordance with the percentage mutated (supplementary_5_). 

.. _supplementary_4: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html
.. _supplementary_5: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html

Extraordinarily, apart from the expected mutation, in all four datasets other pentameers seem to be mutated as well. Especially, CAACC is reported to be highly mutated (around 30%) for all the datasets. This could indicate either the occurrence of other mutations in the dataset or the identification of pentameers which cause error's throughout our workflow.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Variance_occurence_in_MUT.PNG
   :width: 700px
   :height: 450px
   :align: center

   Figure 12: **SNPs analysis for pentameer chunks**. Occurrence of SNP is visualized as a percentage of the amount of time sequence has been reported in the dataset. In all the chunks the middle base has been reported to be mutated in some of the vcf files. The dataset used here is p53 mutated. Interactive figure can be found here(MUT_heatmap_)

.. _script: https://github.com/DouweSpaanderman/NaDA/blob/master/Scripts/vcf_snp_variant_analyser.py
.. _MUT_heatmap: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxMUT_insert_5.0_heatmap_sequences.html

Importantly, the identified alterations could persist anywhere the pentameer is located in the sequence. Furthermore, if an alteration occurs in one specific position it is more prone to be an actual mutation. In contrast, alterations occurring in multiple location with the same pentameer could indicate a systematic problem with pcr or sequencing causing falsely identified mutation. In order to visualize the location of the altered chunks, another script was prepared to compare pentameer locations in the sequence. This script constructed an excel_ with the position of the chunks, including SNP, the percentage of times this particular chunk with SNP occurs in a certain location and the amount of times this specific alteration was found in the dataset. 

As expected, the p53 mutated altered chunk GGATA > T on locations 7578263 to 7578267 was  reported as the only occurrence of this particular sequence with alteration in the whole dataset. Which empowers the conclusion that this is in fact a real mutation. Similarly, other found alterations such as CAACC > G at position 7578333 to 7578337 were also only reported here in the whole sequence. Furthermore, this alteration hasn't been reported yet in GRCh37.13 reference(8). Also, this alteration is located in an intron and is thus less likely to be a driver mutations. Altogether, while it's most likely a passenger mutation due to its consistency between samples, the probability that it's an error cannot be ruled out. 

Currently we are able to identify SNPs occurring in ctDNA when we compare them to cfDNA. Also, we are able to visualize all SNPs occurring in a dataset and are already able to use this knowledge to analyse p53 datasets. However, we have yet to establish enough data, to identify method specific alterations for a broader approach. A future perspective should thus be to gather more cfDNA datasets with Cyclomics, in order to analyse often mutated sequences for data filtering.

Apart from insert data, backbone was also analysed for SNPs. Backbone data could have an application in identifying run specific alterations, which are errors that occur in the sequence and differ in each run. After identification, these run specific errors could be included in data filtering. Firstly, backbone data was similairly stripped and analysed as insert data and also visualized using heatmap. In figure 13 backbone heatmaps from p53 wild-type and 10% mutated were visualized, which clearly shows the occurrence of run specific errors.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/backbone_variance_WT_10.png
   :width: 700px
   :height: 450px
   :align: center

   Figure 13: **SNPs analysis for backbone**. Occurrence of SNP is visualized as a percentage of the amount of time sequence has been reported in backbone. In all the chunks the middle base has been reported to be mutated in some of the vcf files. A) Backbone from p53 wild-type. B) Backbone from p53 10% mutated. Interactive figure can be found here(Backbone_WT_ and Backbone_10_)

.. _Backbone_WT: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_backbone_5.0_heatmap_sequences.html
.. _Backbone_10: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxPool5_xI2_backbone_5.0_heatmap_sequences.html

Similar to insert data, backbone data has also been investigated for SNP locations (excel_). Because backbone shouldn't contain any SNPs, applying reported SNPs to the data filter could provide for run specific background noise canceling. Importantly, SNPs should at least have a 10% occurrence in the backbone to be applied in the data filter, to limit sequence exclusion. 

.. _excel: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Excel.xlsx

Script Tests
++++++++++++
Before functions and scripts are run over multiple files and directories, they should be checked for quality. In order to check a function for its functionality, test scripts can be written. These testing scripts use the assert function to identify if the set criteria are met. As an example the earlier described parse_fasta_file_error function is checked for its quality with the following testing script::

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

The class function is used to define which function is going to be checked for quality. Firstly the function is setup with a test file, this file is designed to identify flaws in the function. In other words, it consists off alot of errors which the function should not pickup. Next, multiple assertions are made, such as the assertion that letters in sequence can only consist of A, C, T, G and N. Also score should consist of characters and not involve any letters. While this is an example of a test script, every function and script is checked for there quality. Testing scripts can be found here_. 

.. _here: https://github.com/DouweSpaanderman/NaDA/tree/master/Testing%20scripts

|
|

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Next.png
   :align: center
   :width: 100px
   :height: 100px
   :target: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Conclusion.html