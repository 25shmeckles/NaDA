Introduction
------------
Tumor genotyping is essential for improving the clinical outcome for many cancer patients (1_). Not only it is important for tailoring drugs specific targets for the genetic landscape of the tumor, it is also vital for better diagnostics and for follow up monitoring. Currently, tumor genotyping is conducted using samples of local biopsies of tumor tissue. Nonetheless, this method of identifying the genetic landscape of malignancies is clearly limited by the heterogeneity present in the samples, that exists as a result of tumor infiltration in surrounding tissue and central necrosis due to rapid tumor growth. Furthermore, tissue biopsies can only be conducted a limited number of times, which poses a problem for accurate identification of the tumors genetic landscape (2_). Therefore, new minimal invasive methods are required to be able to accurately identify the tumors genetic landscape.

Just recently, circulating tumor DNA (ctDNA) has been suggested as a new innovating method for identifying the tumor genetic landscape(3_). These DNA fragments are derived from necrotic tumor cells and subsequently further fragmented by macrophages (4_). ctDNA makes up for a small percentage of the whole cell-free DNA (cfDNA), increasingly so upon heavy tumor burden. Theoretically, ctDNA provides for a better alternative for tumor genotyping, because of its minimal invasive nature and therefore ability to be used frequently. Most important of all it is not limited by tumor heterogeneity, as it accounts for both intratumoral and intermetastatic and has been shown to contain multiple types of alterations. However, the use of ctDNA comes with its own set of difficulties. First, the low quantity allelic fractions of ctDNA in patient's blood(5_), especially in patients with a low tumor burden, as it is simply harder to yield sufficient ctDNA for performing genetic tests. Also, discrimination between cfDNA and ctDNA remains problematic. These issues have to be addressed and sufficient improvements have to be made by researchers to develop highly specific and sensitive clinical tests in the future.

In order to tackle these hurdles, The Kloosterman group (UMC Utrecht) uses a new and innovating method called Cyclomics to identify alterations in ctDNA (Figure 1). Current research on this method is conducted to detect p53 mutations in patients with head and neck carcinomas. When a sensitive and specific test for p53 mutations has been devised, Cyclomics could be applied in the identification of mutations in other driver mutations. Cyclomics uses rolling circle amplification to increase low allelic fractions of mutated ctDNA fragments(6_). This method of enhancing the amount of ctDNA is known as library preparation, which causes Nanopore sequencing to increases in sensitivity. However, as with all sequencing methods, sequencing is not flawless and covers real biological mutations, but also asymmetric DNA errors, PCR and sequencing errors (7_). Therefore, including multiple filtering steps upon the data achieved from the Nanopore sequencer is essential for identifying real biological mutations and thus improving specificity and sensitivity of the clinical tests. In this report, data analyzing will be covered with Python, identifying regions with high and low mutational value as well as identifying single nucleotide polymorphisms (SNPs) in a p53 dataset. Our analysis is able to identify real driver mutations by comparing wild-type and mutant cfDNA on p53. Also, wild-type analysis shows the occurrence of both real passenger mutations and method induced errors, which can be used in future data filter construction.
 
.. figure::  https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Figure_workflow.jpg
   :width: 840px
   :height: 450px

   Figure 1: **Workflow for Cyclomics**. A) for library preparation, patient blood is isolated, which contains circulating tumor DNA (ctDNA) fragments. Using a phosphorylated backbone, these fragments can be circularized and subsequently amplified with rolling circle amplification. This leaves very long reads, which are ideal for nanopore sequencing. B) Nanopore sequencing is used to identify nucleotides in a read and after base and variant calling gives fastq and vcf files. After applying a Data filter, true driver mutation can be identified. Alterations visible in the ctDNA can be separated by real driver mutations, Real passenger mutations, asymmetric DNA errors and PCR and sequencing errors. Driver mutations are alterations in a gene that gives a selective growth advantage, while passenger mutations do not. Errors can be caused by polymerase induced errors in PCR or wrong base calling in sequencing. Finally, asymmetric DNA errors are alteration which is only present in one strand.

.. _1: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/References.html
.. _2: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/References.html
.. _3: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/References.html
.. _4: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/References.html
.. _5: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/References.html
.. _6: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/References.html
.. _7: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/References.html

|
|

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Next.png
   :align:  center
   :width: 100px
   :height: 100px
   :target: https://rawgit.com/DouweSpaanderman/NaDA/master/Documentation/build/html/General%20concept.html