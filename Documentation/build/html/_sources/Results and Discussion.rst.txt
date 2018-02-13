Results and discussion
----------------------
Fastq files
+++++++++++
**Quality score distribution**. Fastq is a text-based format derived from sequencers such as Illumina and Nanopore that contain a biological sequence and for every base a quality score.
this quality score represents the sequencers ability to call a base as true, in other words, the score gives an indication about how sure the sequencer is that the output it gives is the base that is actually present in the sample. If the score for the base is low, it's most likely the right base, while if the score is high the machine is unsure if it's the real base. Nanopore sequencer measures changes in current to identify bases in the sample and if changes in this current are less clear a high quality score is given for this base.

Analyzing quality scores on sequences is important to be able to understand which regions of bases are difficult to get base called. More importantly, a mutation called in a region with low quality score is more likely to be a real mutation in contrast to a mutation that occurs through wrong base calling in areas with high quality scores. Clarifying which regions have high scores is important for achieving tests with high specificity, as these regions will not be selected for variant calling. Likewise, regions with a low score are good targets for identifying mutations and are thus able to increase test sensitivity. Therefore the aim of fastq file analysis is to identify regions of bases with high and low scores, as these regions can be used for data filtering. In data filtering, high score regions can be excluded and low score regions can be included because mutations identified in high score regions are more likely to be sequencing errors. 
 
In order to analyze high and low score regions, fastq files first need to be stripped of bases and their scores. the following function was written to strip the files of data:: 

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

The function returns a dictionary with identification (id), sequence and scores for individual bases in the fastq file. With the data stripped from the fastq file, the first question was whether there are changes in quality score present in our file. Therefore, the distribution of the error rate was visualized for every individual base in the sequence as shown in Figure 2, which illustrates the occurrence of high and low score regions. In these high score regions, bases have more chance to be falsely assigned and are thus less reliable for mutation identification. Likewise, mutations found in low score regions are more likely to be rightly assigned and thus true mutations. Figure 2 clearly visualizes changes in quality score distribution and could thus indicate that high and low regions are present in our sample.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Fastq_files_qualityscore.png
   :width: 840px
   :height: 560px

   Figure 2: **Distribution of error rate (quality score) for every individual base.** In this Figure, every 10th base is visualized on the X-axis to clarify the data. The mean of all the quality scores is 0.21.

While Figure 2 gives a good overview on the importance of analyzing quality scores, one fastq file has little quantifiable value for identifying high score regions. To further investigate our first question, whether changes in quality scores are present in our samples, 1139 fastq files from Nanopore with reads on chromosome 9 were likewise investigated for quality score distribution. Firstly, mean scores of the quality scores were plotted in Figure 3. The mean scores of the files are a good first indication of the sequence quality and can be used to filter and select files for variant calling, in which a low mean quality score is important. Similar to Figure 2, mean score distribution also visualizes changes in quality scores between files further empowering the possibility that high and low regions are present in our samples.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Mean_distribution_of_all_quality_score.png
   :width: 840px
   :height: 820px

   Figure 3: **Mean distribution of quality scores**. Mean scores from 1139 fastq files were calculated and visualised in a violin plot.

**Region quality score distribution**. Now we know that quality score changes occur in our samples, the next question arose whether these changes translate into regions with high and low scores. In order to investigate if these regions are present in our sample, regions of bases were selected instead of single bases. Also, these regions were investigated in different sizes to analyze the impact of different sized regions on the occurrence of high and low score region. Therefore, sequences were chunked into pieces of four (tetramer) to six (hexamer) bases and for every sequence chunk, its respected mean score was calculated. Size and overlap of the chunks of sequences could be selected by the following function::

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

Following, the mean of all the quality scores of the same chunks of sequences were either plotted directly (Figure 4 A - C) or indirectly after being divided into categories of high, medium and low quality score (Figure 4 D - F). Categorizing was performed before calculating the mean of a sequence, aiming to manipulate and increase data analysis. Subsequently, categories were counted for each sequence. Parameters for categorizing were randomly selected and differentiated for each size, because with larger regions the mean of the quality score shift further towards medium. To account for this difference the high requirements are lowered and the low requirements were increased as follows::

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
   :width: 840px
   :height: 554px

   Figure 4: **Quality score analysis in 6 scenarios.** A - C) Mean score for all combinations in chunk sizes (A = 4, B = 5, C = 6) for 1139 fastq files derived from nanopore sequencing of chromosome 9. D - F) Scores for regions have been categorized into high, medium and low for regions of the same size as A to C. Next, the number of times a region was called under a certain category was counted and collected for the same data set. In these figures, scores are set in percentage of the total amount of times a region occurs in the data set.(Interactive figure at GridPlot_)

.. _GridPlot: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/gridplot.html

Table 1 highlights the highest and lowest five scoring sequences measured in Figure 4. The highest scoring sequence has the biggest chance to have wrongly assigned bases in its sequence. In contrast, bases in lower scoring sequences have been most likely rightfully assigned and are therefore indeed the right base. These findings should be taken into account when investigating mutations, as a mutation found in, for instance, TTCC is more likely to be a real mutation than an alteration found in GCTT. Conclusively, this dataset of chromosome 9 does indeed have high and low quality score regions which should be further analyzed for respectively exclusion and inclusion into a data filter. 

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

**Clustering**. Overall, the fastq analysis showed the occurrence of high and low regions in our dataset. Next, these regions should be further investigated and clearer visualized for construction of our data filter. In order to analyze which regions should be included and excluded in our data filter clustering is conducted. Clustering is a method in which data points get coupled in groups (clusters) by a certain geometry. Here, K-Means is used for clustering, which makes clusters based on the distance between data points and a selected centroid. The centroid is the mean of a cluster and is defined by a trial and error process. This process is repeated until centroids are selected, which happens when the within-cluster sum of squares is minimized.

In Figure 5 three clusters are formed in which the yellow one represents sequences with often reported high score and few times reported low score. In this cluster, alterations are more likely to be falsely assigned. Furthermore, in the blue cluster, alterations are more likely to be rightfully assigned, cause sequences are often reported in low score and only a few times in high score. A great advantage of clustering of data is that a clearer view on which sequences to include and exclude for mutation calling can be provided. Here, regions in the yellow cluster should be excluded for variant calling. In other words, any alteration found in regions located in the yellow cluster should be filtered. 

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/clusterplot.png
   :width: 840px
   :height: 325px

   Figure 5: **Clustering of data from 1139 fastq files.** Sequences of 5 bases (pentamer) are measured for quality score and reported in high or low score. Here, the percentage of the number a sequence is reported in high and low score is visualized. Following, clustering was conducted using K-means, thereby producing three clusters. 

Together, fastq data indicates that quality score is important for identifying regions which are promising for mutation calling and which regions should be avoided. As described earlier, regions with a high quality score should be avoided while looking for mutations. In contrast, low quality score regions have the potential for identifying real mutations in ctDNA. For the construction of a data filter, regions should be selected from the clustering data. Nevertheless, only fastq files of one run on chromosome 9 have now been investigated and while this can provide for data filtering of this specific run, for data filter construction on a larger scale, multiple runs should be compared. Importantly, the quality scores of sequences can differ on the method being used. In this case, our method involves rolling circle amplification and nanopore sequencing of cfDNA. In order to make a sensitive data filter, a big database of healthy cfDNA should be investigated on quality score for sequences. In this way, the filter can exclude and include regions with respectively high and low quality scores. Furthermore, quality scores of sequences can differ on every run, causing some sequences to have higher or lower scores. This discrepancy should also be accounted for in the data filter. A possibility is to either include healthy cfDNA into every run or compare backbone sequences to identify run specific sequence quality score differences. Backbone analysis for identifying run specific changes in quality score sound more promising as the backbone is already incorporated in sequencing and the backbone sequence is known. 

**p53 wild-type and mutant dataset analysis**. So far, only run specific sequence quality scores have been discussed. However, our goal is to identify method specific high and low quality score regions. In order to visualize high and low score sequences specific for our method, data analysis should be conducted on multiple runs. Analysis of multiple runs can be simultaneously conducted using the High-Performance Computing (HPC) facility in the UMC, which will be discussed in the following segment. For these multiple datasets, aims are similar to the chromosome 9 dataset. Firstly, the mean should give an indication on whether changes in quality score occur in these datasets. Secondly, score plotting could tell us if regions exist. Finally, clustering would provide us with groups of regions to include and exclude in the data filter construction.

On the HPC, multiple ctDNA datasets derived from cyclomics are available for analysis, as they have been earlier constructed by the Kloosterman group. Here, the focus goes towards the rolling circle amplification of p53 mutated and wild-type(WT) datasets. The p53 gene codes for a tumor suppressor protein that can initiate apoptosis, arrest in the cell cycle and activate DNA repair proteins. In order to analyze these p53 datasets, a fastq_script_ was written to achieve similar data analysis and visualization as was applied to fastq files from chromosome 9. To improve data visibility, minor visualization updates were conducted. 
Both mutated and WT datasets were separated in equally sized chunked files, around 4000 files each, and analysed as individual chunks to increase script parallelization, thus increasing speed. First, the datasets were investigated for changes in quality score, to give an indication whether these changes also harnessed in these datasets similarly to the chromosome 9 dataset. Therefore, the mean scores of all files were calculated and visualized in Figure 6.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/RCAxMUT_WT_boxplot.png
   :width: 840px
   :height: 840px

   Figure 6: **Boxplot of mean score from p53 mutated and WT dataset.** For both datasets chunk 0 to 9 have been visualized. Chunk 10 to 24 were excluded, but showed similar results.

.. _fastq_script: https://github.com/DouweSpaanderman/NaDA/blob/master/Scripts/fastq_qualityscore_analyser.py

This boxplot clearly visualizes the lack of consistency between quality scores in the same sequence run. Similarly to our previous chromosome 9 dataset, Figure 6 visualizes the occurrence of difference in quality score in p53 datasets. However, the inconsistency between the different chunks in the same dataset reveals that quality scores are highly fluctuating between chunks. This could indicate that regions reported in high and low quality score also shift heavily between chunks. 

While mean scores give an indication on quality score analysis, both quality score plotting and clustering are yet to determine whether high and low score regions exist and persist in multiple chunks and datasets. For every chunk derived from a dataset, sequences have been analyzed and visualized as tetramer, pentamer and hexamer. This was conducted to analyze the impact of different sized regions on the occurrence of high and low score region. Here, tetramers of wild-type chunk 0 to 3 have been plotted as shown in Figure 7.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Fastq_gridplot_WT.png
   :width: 840px
   :height: 680px

   Figure 7: **Tetramer sequence analysis for chunks zero to three of the p53 wild-type database.** This figure illustrates the mean quality score for each tetramer possible in one chunk. A) chunk 0. B) chunk 1. C) chunk 2. D) chunk 3. The interactive figure can be found here and also visualizes data analysis when divided into high, medium and low group.(WT_chunk0_, WT_chunk1_, WT_chunk2_ and WT_chunk3_)

.. _WT_chunk0: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk0_4.0_3.0_score_plotting.html
.. _WT_chunk1: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk1_4.0_3.0_score_plotting.html
.. _WT_chunk2: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk2_4.0_3.0_score_plotting.html
.. _WT_chunk3: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk3_4.0_3.0_score_plotting.html

Similarly to the boxplot, there seems to be a lack of consistency between chunks as high reported tetramers differ heavily between these chunks. Similar results are visible for bigger sized chunks(supplementary_1_) and chunks derived from the p53 mutant dataset(supplementary_2_). These datasets show that there is yet to be proven that a correlation between quality scores and specific regions or chunks exists. However, clustering could clarify the occurrence of high quality score regions by better identification of these regions. In order to cluster data derived from dataset chunks, the same algorithm (K-Means) is used as described above. In Figure 8 clustering of chunks, 0 to 3 from p53 WT has been visualized.

.. _supplementary_1: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html
.. _supplementary_2: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Fastq_gridplot_WT_cluster.png
   :width: 840px
   :height: 815px

   Figure 8: **Clustering of hexameer sequence for chunks zero to three of the p53 wild-type database.** Data points are visualised as percentage reported in high (y-axis) and low score(x-axis). A) chunk 0. B) chunk 1. C) chunk 2. D) chunk 3. Interactive cluster plot can be found here. (WT_chunk0_cluster_, WT_chunk1_cluster_, WT_chunk2_cluster_ and WT_chunk3_cluster_)

.. _WT_chunk0_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk0_6.0_5.0_score_clustering.html
.. _WT_chunk1_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk1_6.0_5.0_score_clustering.html
.. _WT_chunk2_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk2_6.0_5.0_score_clustering.html
.. _WT_chunk3_cluster: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_chunk3_6.0_5.0_score_clustering.html

By clustering, regions can be identified that have both been reported often as high score and few times as low scores. However, between chunks there is a huge discrepancy in quality scores. Chunk 0 and 2 have an overall much lower quality score in comparison with chunk 1 and 2, which was also identified with mean quality scores in Figure 6. This big difference in overall quality of the chunks suggests that quality score cannot be used for constructing a data filter as regions vary too much. Nevertheless, if regions are always present in the same cluster between chunks, high and low quality score regions could still be identified. Additionally, quality score means should not have to influence score from a single region in comparison to other regions in the same chunk. With this in mind, future analysis should be conducted comparing regions between chunks for the same group. For example, regions reported in the high score group from chunk 0 should be compared to regions reported in the high score groups of the other chunks.

As described earlier, identifying high quality score regions could help data filtering for mutations, as bases in regions with high quality score are less likely to be rightfully assigned and the other way around for low quality score region. Therefore, respectively, these regions could be excluded or included in variant calling. Currently, quality score analysis shows a lot of inconsistency between chunks of the same dataset. Thus, it is yet impossible to conclude any regions that have either a high or low quality score. Therefore, quality score has currently no application in creating a data filter. Altogether, although quality score shows limited possibilities, regions should still be compared between chunks, because regions could still be reported in similar clusters (high score or low score clusters). Further comparison of regions between chunks could clarify whether there is a correlation between regions and quality scores. 

Variant Call Format files
+++++++++++++++++++++++++
**Mutation distribution of single nucleotide polymorphisms**. After base calling is conducted, reads are compared with the reference genome in order to call variants occurring in the samples. The output of variant calling is located in Variant Call Format (VCF) files, which are in text-based formats similar to fastq, containing data of single positions bases in the genome. In these files, variants are formatted with the reference included. For sequenced sites, the number of reads found with the mutation and reference is given.

In these VCF files multiple variances can be identified depending on which dataset is used. These variances identified in the VCF files can be real mutations but also PCR and sequencing errors. Analysing VCF files could yield information about real mutations positions as well as the positions where errors are more prone to be located. Similarly to fastq files, it could be possible that errors occur more often in specific regions of bases. Therefore, Alterations found in these locations could be excluded by a data filter. Our goals for VCF analysis are to be able to identify real driver mutations and to analyze if regions can be identified which are more prone to errors. Identifying real driver mutations can be done by comparing wild-type (WT) and mutated (Mut) databases as driver mutations should not be included in the WT. Furthermore, region analysis can be conducted by constructing larger sequences similar to fastq analysis. 

In order to achieve our goals, a dataset was analyzed derived from the Cyclomics project, in which sequencing was performed with Nanopore. The data contains a part of the p53 gene on chromosome 17, existing of 160 nucleotides and a backbone, which was used for the circular PCR reaction. In total, 1187 VCF files were used for variant calling. These VCF files were screened for single nucleotide polymorphism (SNP) occurrence. Files were stripped of reported mutated bases and other data was discarded. As described earlier, every variant site has a number of reads that cover this site. These reads can be both coupled to the mutation and the reference. For example, on position 7577503 a SNP was found in 6 reads and 3 reads were coupled to the reference. While normally the number of reads coupled to the mutation in contrast to the total amount of reads is important, here the occurrence of certain SNPs has been firstly investigated. Analysing SNPs in ctDNA could help identify real driver mutations present in a tumor. Likewise, analyzing cfDNA for SNPs by identifying passenger mutation, asymmetric DNA errors, PCR errors and sequencing errors can be useful for lowering background noise.

Our first aim is to analyze and visualize SNP occurrences in the p53 databases. Identifying often manifesting alteration can help understand what bases are more prone to errors. Normally, no discrepancy between alterations should occur as biological mechanisms do not favor specific base substitutions. Therefore, changes in Cyclomics can be devised to suppress these errors. In order to investigate the amount of SNPs in the files, VCF files were similarly stripped as Fastq files and separated by either sequence or backbone. Next, for the alterations, a parameter was set at a minimum of 25 percent of the reads that should be coupled to the mutant variant and visualized in Figure 9:

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Combined_vcf_snp_analysis.png
   :width: 840px
   :height: 405px

   Figure 9: **Distribution of SNPs in the sequence of the p53 gene for 1187 VCF files.** Parameter for variant identification was set at 25% of the reads to the variant. Variants are displayed as C > T, meaning that T substitutes C. A) Bar plot with single nucleotide polymorphisms occurrence as percentage of whole. B) Heatmap of the same variances with the number of occurrences in the files

Both Figure 9A and 9B illustrate a common occurrence of G > A mutation and to lesser extend C > A. The high prevalence of these SNPs in contrast to other alterations are a strong indication that these alterations are caused by a non-biological mechanism, which can be errors in the rolling circle amplification, library preparation and sequencing of the ctDNA. In literature, cytosine deamination has been described to increase C:G > T:A noise levels (8_). Cytosine deamination occurs when the amino group of cytosine is removed resulting in a change to a uracil analog. Also, the less occurring alteration C > A has been reported to be caused by oxidative DNA damage during sample preparation(9_). Both these types of alterations can be a result of polymerase-induced errors. Possible suggested methods to suppress these errors are adding DNA repair mechanisms upon polymerase chain reaction (PCR) and lowering heat. However, if results are similar in other databases, an in silico approach to polish background noise can also be devised. Overall, optimizing Cyclomics to minimize error occurrences should first be conducted.

.. _8: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _9: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html

**Region SNP analysis**. As described earlier, one of our goals is to identify if regions are present in our sample which are more prone to errors. Therefore, trimers and pentamers were selected in which the middle base was reported to contain a SNPs in some of the files for heatmap analysis. Creating chunks for SNP analysis instead of single base analysis was conducted to visualize sequences that were more likely to contain SNPs. Pandas library was used to create a dataframe for the number of times mutation occurred to either A, T, C or G. This dataframe was then mapped to a heatmap with the reference sequence (Figure 10). Just as in previous Figures, length of the surrounding bases can be changed to give a wider variety of information. Firstly, analyzing of this p53 dataset visualized the occurrence of regions that are prone to errors. Secondly, This Figure gave more information about base combinations with high alteration affinity, such as CGC to CAC.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Variance_occurence_in_sequence_vcf_3.PNG
   :width: 840px
   :height: 335px

   Figure 10: **Occurence of variance per reference sequence to different bases.** In all the sequences the middle base is reported to be mutated in some of the vcf files. This mutation again has a parameter that is set at 25% of the reads at least mutated. 

Identifying high variance chunks in healthy cfDNA is important to understand where errors arise in the sequence and thus could help polishing background noise by deselecting these chunks for variant calling analysis. However, it is important to keep in mind that variance present in WT samples could not only be errors but also real passenger mutations. 

Furthermore, just as with the fastq files, variances can be separated between alterations specific for a run and alterations specific for the method being used. For instance, CTC > A could be an alteration that is specifically highly mutated in a particular run, while CGC > A occurs often in every run with this method of rolling circle amplification and Nanopore sequencing. Therefore, filtering should be able to account for both run specific and method specific alterations. In the same manner, a large database of healthy cfDNA could accomplish a method specific filter and adding healthy cfDNA into every run could provide with a specific alteration filter, which will be helpful in lowering background noise. Also more convenient, backbone data could be used to identify run specific errors as the backbone does not change between runs and should thus never contain alterations.

**p53 wild-type and mutant dataset analysis**. In the previous segment, only one run was investigated for the presence of SNPs. In order to establish enough information for specifying a method specific data filter, multiple Cyclomics databases should be analyzed for occurring alterations. In the following segments analysis of four datasets was conducted. These datasets all cover the p53 gene, either WT, Mut, 1% mutated or 10% mutated. Investigation of these datasets could provide for real driver mutations and often error harnessing regions. Further analysis on these VCF files was conducted on the HPC system as it required greater computer performance and parallelization. In order to conduct VCF analysis on the HPC, a VCF analysis script_ was written. This script analyses mutation occurrences in the insert and backbone similarly as described above. Subsequently, visualization was updated, presenting sequences as a percentage of times it has been reported to contain a SNP to how many times it was reported in total for the whole dataset. 

Firstly, SNPs were investigated for their individual positions. This was conducted to be able to identify the real driver mutation by comparing wild-type and mutated datasets. As visualized in Figure 11, using this single position heatmap, a clear mutation, position 7578265 A > T, can be found in the p53 mutated dataset. 

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/single_base_insert_MUT.PNG
   :width: 840px
   :height: 325px

   Figure 11: **Single position SNPs analysis for p53 mutated database**. Occurrence of SNP is visualized as a percentage of times it has been reported with a SNP to the total amount of times it has been reported in the dataset. Interactive figure can be found here(single_base_insert_MUT_)

.. _single_base_insert_MUT: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxMUT_single_base_insert_1_heatmap_sequences.html

Similar to Figure 11, wild-type, 1% and 10% have been analyzed for position specific SNPs (supplementary_3_). As expected, wild-type shows no occurrence of the specific mutation and the other datasets are in concordance with the percentages of reads that should contain the SNP. Therefore, this script_ is able to identify real driver mutations in this p53 dataset, which is A to T at position 7578265. 

.. _supplementary_3: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html

However, while the expected alteration is found in our datasets, for those datasets other alteration were also identified. These alterations were not in line with our expectations and could have a number of explanations. Firstly, it could be a real passenger mutation which is not normally present in the reference. Secondly, these alterations could be caused by protocol errors, such as PCR errors and sequencing errors. For both these types of alterations could be accounted for by respectively excluding certain positions and certain regions of bases in the data filter. However, the distinction between these alterations is difficult to determine. Currently, found alterations can help to construct a data filter specific for p53 databases as specific positions can be excluded. Nevertheless, our aim is to construct a data filter that can be applied to multiple gene databases. Therefore, regions should further be investigated to determine if they are reported as mutated because of errors or real passenger mutations. 

In order to construct a data filter usable for multiple genes, regions rather than single positions were investigated for SNPs. Region selection could help identifying SNPs due to protocol errors and help exclude them from the analysis. First, chunks of trimers and pentamers were selected as described in region SNP analysis. Subsequently, these chunks were plotted in a heatmap visualized in Figure 12. Mutation occurrence (GGATA > T) could again be clearly visualized in this dataset. Furthermore, WT analysis showed a complete absence of this mutation (supplementary_4_) and the other two datasets are in concordance with the percentage mutated (supplementary_5_). 

.. _supplementary_4: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html
.. _supplementary_5: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Supplementary.html

Extraordinarily, apart from the expected mutation, in all four datasets other pentamers seem to be mutated as well. Especially, CAACC is reported to be highly mutated (around 30%) for all the datasets. This could indicate either the occurrence of other mutations in the dataset or the identification of pentamers which cause errors throughout our workflow.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Variance_occurence_in_MUT.PNG
   :width: 840px
   :height: 325px

   Figure 12: **SNPs analysis for pentameer chunks**. Occurrence of SNP is visualized as a percentage of the amount of time sequence has been reported in the dataset. In all the chunks the middle base has been reported to be mutated in some of the vcf files. The dataset used here is the p53 mutated dataset. The interactive figure can be found here(MUT_heatmap_).

.. _script: https://github.com/DouweSpaanderman/NaDA/blob/master/Scripts/vcf_snp_variant_analyser.py
.. _MUT_heatmap: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxMUT_insert_5.0_heatmap_sequences.html

Importantly, the identified alterations could be present anywhere the pentamer is located in the sequence. Furthermore, if an alteration occurs in one specific position it is more prone to be an actual mutation. In contrast, alterations occurring in multiple locations within the same pentamer could indicate a systematic problem with PCR or sequencing causing a falsely identified mutation. In order to visualize the location of the altered chunks, another script was prepared to compare the pentamer locations in the sequence. This script constructed an excel_ with the position of the chunks, including SNP, the percentage of times this particular chunk with SNP occurs in a certain location and the number of times this specific alteration was found in the dataset. 

As expected, the p53 mutated altered chunk GGATA > T on locations 7578263 to 7578267 was reported as the only occurrence of this particular sequence with alteration in the whole dataset. This finding empowers the conclusion that this is, in fact, a real mutation. Similarly, other found alterations such as CAACC > G at position 7578333 to 7578337 were also only reported here in the whole sequence, suggesting that this is also a real mutation. However, this alteration has not been reported yet in the GRCh37.13 reference genome(10_), questioning the importance of this mutation. Also, this alteration is located in an intron and is thus less likely to be a driver mutation. Altogether, a mutation like CAACC > G is most likely a real passenger mutation and this location should thus be excluded by our data filter. Besides this passenger mutation, our datasets also show signs of error-prone regions. In figure 12, alteration CTGGG > A was found with an occurrence of 7,5%. Importantly, location analysis described in excel_ showed the presence of this alteration in multiple locations as both position 7578321 and 7578330 harnessed this specific alteration. Therefore, CTGGG > A is a region that is most likely more prone to errors occurring in our method and should be excluded by our data filter. 

Future investigations should be conducted towards automation of this process of identifying reported alterations by error-prone regions or real passenger mutations. Furthermore, automation of this process should also check if the mutation is reported in the GRCh37.13 reference genome(10_). Automation of this process is important for further investigation of VCF files for data filter construction.

.. _10: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html

Currently, we are able to identify SNPs occurring in ctDNA when we compare them to cfDNA. Additionally, we are able to visualize all SNPs occurring in a dataset and are already able to use this knowledge to analyze p53 datasets. Also, error-prone regions have been shown to persist in our database and can be further investigated for data filter construction. As a future perspective identification of both passenger mutation and error-prone regions should be automated to make data filter construction easier.

Apart from insert data, the backbone was also analyzed for SNPs. Backbone data could be applied in identifying run specific alterations, which are errors that occur in the sequence and differ in each run. After identification, these run specific errors could be included in data filtering. Firstly, backbone data was similarly stripped and analyzed as insert data and also visualized using heatmap. In Figure 13 backbone heatmaps from p53 wild-type and 10% mutated were visualized, which clearly shows the occurrence of run specific errors.

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/backbone_variance_WT_10.png
   :width: 840px
   :height: 665px

   Figure 13: **SNPs analysis for backbone**. Occurrence of SNP is visualized as a percentage of the amount of time sequence has been reported in backbone. In all the chunks the middle base has been reported to be mutated in some of the vcf files. A) Backbone from p53 wild-type. B) Backbone from p53 10% mutated. Interactive figure can be found here(Backbone_WT_ and Backbone_10_).

.. _Backbone_WT: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxWT_backbone_5.0_heatmap_sequences.html
.. _Backbone_10: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/source/_static/RCAxPool5_xI2_backbone_5.0_heatmap_sequences.html

Similar to the insert data, the backbone data has been investigated for SNP locations (excel_). Because the backbone should not contain any SNPs, applying reported SNPs to the data filter could provide for run specific background noise canceling. Importantly, SNPs should at least have a 10% occurrence in the backbone to be applied in the data filter, to limit sequence exclusion. Furthermore, automation should be conducted to analyze error-prone regions.

.. _excel: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Excel.xlsx

Script Tests
++++++++++++
Before functions and scripts are run over multiple files and directories, they should be checked for quality. In order to check a function for its functionality, test scripts can be written. These testing scripts use the assert function to identify whether the set criteria are met. As an example the earlier described parse_fasta_file_error function is checked for its quality with the following testing script::

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

The class function is used to define which function is going to be checked for quality. Firstly, the function is setup with a test file, which is designed to identify flaws in the function. In other words, the test file consists of a lot of errors which the function should not pickup. Next, multiple assertions are made, such as the assertion that letters in sequence can only consist of A, C, T, G and N. Also, the score should consist of characters and not involve any letters. While this is an example of a test script, every function and script is checked for their quality. Testing scripts can be found here_. 

.. _here: https://github.com/DouweSpaanderman/NaDA/tree/master/Testing%20scripts

|
|

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Next.png
   :align: center
   :width: 100px
   :height: 100px
   :target: https://rawgit.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Conclusion.html