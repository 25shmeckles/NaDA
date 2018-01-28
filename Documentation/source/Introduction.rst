Introduction
------------
Tumor genotyping is essential for improving the clinical outcome for many cancer patients. Not only is it important for tailoring drugs specific targets for the genetic landscape of the tumor. It is also vital for better diagnostics as well as follow up monitoring. Furthermore, increasing knowledge of tumor tissue genotyping can enhance our understanding of the molecular mechanisms of cancer. Today, tumor genotyping is conducted using local biopsies of tumor tissue. This method of achieving the genetic landscape is, however heavily limited mostly due to its invasive nature. Therefore, tissue biopsies can only be conducted a limited amount of times, enhancing difficulties with tumor heterogeneity and the fluctuation of the genetic landscape of the tumor in time (1_). Altogether, these limitations require a new none invasive method to be able to tackle these difficulties.

Just recently, circulating tumor DNA (ctDNA) has been suggested as a new innovating method for identifying the tumor genetic landscape. These DNA fragments are derived from necrotic tumor cells and subsequently further fragmented by macrophages (2_). ctDNA makes up for a portion of the whole cell-free DNA (cfDNA), increasingly so upon heavy tumor burden. Theoretically, ctDNA provides a better alternative for tumor genotyping, because of its noninvasive nature, ability to be used unlimited times, it accounts for most of both intratumoral and intermetastatic heterogeneity and has been shown to contain multiple types of alterations. However, ctDNA encounters some difficulties, such as low allelic fractions of ctDNA in patient's blood(3_), especially on low tumor burden. Also, discrimination between cfDNA and ctDNA has seemed to be problematic. These problems hinder researchers to develop high specific and sensitive clinical tests.

The Kloosterman group (UMC Utrecht) uses a new and innovating method called Cyclomics to identify alterations in ctDNA (Figure 1). Current research upon this method is conducted for p53 mutation in patients with head and neck carcinoma's, but upon achieving a sensitive and specific test, Cyclomics could provide for identifying multiple driver mutations in ctDNA. Cyclomics uses rolling circle amplification to increase low allelic fractions of mutated ctDNA fragments(4_). This library preparation for Nanopore sequencing increases its sensitivity. However, as with all sequencing methods, sequencing isn't flawless and covers real biological mutations, but also asymmetric DNA errors, PCR and sequencing errors (5_). Therefore, including multiple filtering steps upon the data achieved from the Nanopore sequencer is essential for identifying real biological mutations and thus improving specificity and sensitivity of the clinical tests. In this report, data analyzing will be covered with Python, identifying regions with high and low mutational value as well as identifying single nucleotide polymorphisms (SNPs) in a p53 dataset.
 
.. figure::  https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Figure_workflow.jpg
   :width: 600px
   :height: 400px
   :align: center

   Figure 1: **Workflow for Cyclomics**. A) for library preparation, patient blood is isolated, which contains circulating tumor DNA (ctDNA) fragments. Using a phosphorylated backbone, these fragments can be circularized and subsequently amplified with rolling circle amplification. This leaves very long reads, which are ideal for nanopore sequencing. B) Nanopore sequencing is used to identify nucleotides in a read and after base and variant calling gives fastq and vcf files. After applying a Data filter, true driver mutation can be identified.

.. _1: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _2: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _3: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _4: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _5: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html

|
|

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Next.png
   :align: center
   :width: 200px
   :height: 100px
   :target: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/General%20workflow.html