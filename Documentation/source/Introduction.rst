Introduction
------------
Tumor genotyping is essential for improving the clinical outcome for many cancer patients (1_). Not only it is important for tailoring drugs specific targets for the genetic landscape of the tumor. It is also vital for better diagnostics and for follow up monitoring. Today, tumor genotyping is conducted using local biopsies of tumor tissue. This method of achieving the genetic landscape is, however heavily limited mostly due to its invasive nature. Tissue biopsies can only be conducted a limited amount of times, causing difficulties with tumor heterogeneity and the fluctuation of the genetic landscape of the tumor in time (2_). Altogether, these limitations require a new minimal invasive method to be able to tackle these difficulties.

Just recently, circulating tumor DNA (ctDNA) has been suggested as a new innovating method for identifying the tumor genetic landscape(3_). These DNA fragments are derived from necrotic tumor cells and subsequently further fragmented by macrophages (4_). ctDNA makes up for a portion of the whole cell-free DNA (cfDNA), increasingly so upon heavy tumor burden. Theoretically, ctDNA provides a better alternative for tumor genotyping, because of its minimal invasive nature, ability to be used often, it accounts for most of both intratumoral and intermetastatic heterogeneity and has been shown to contain multiple types of alterations. However, ctDNA runs into some difficulties, such as low allelic fractions of ctDNA in patient's blood(5_), especially on low tumor burden. Also, discrimination between cfDNA and ctDNA has seemed to be problematic. These problems hinder researchers to develop high specific and sensitive clinical tests.

The Kloosterman group (UMC Utrecht) uses a new and innovating method called Cyclomics to identify alterations in ctDNA (Figure 1). Current research on this method is conducted to detect p53 mutations in patients with head and neck carcinomas. When a sensitive and specific test for p53 mutations has been devised, Cyclomics could be applied in the identification of other driver mutations. Cyclomics uses rolling circle amplification to increase low allelic fractions of mutated ctDNA fragments(6_). This method of enhancing the amount of ctDNA is known as library preparation, which causes Nanopore sequencing to increases in sensitivity. However, as with all sequencing methods, sequencing isn't flawless and covers real biological mutations, but also asymmetric DNA errors, PCR and sequencing errors (7_). Therefore, including multiple filtering steps upon the data achieved from the Nanopore sequencer is essential for identifying real biological mutations and thus improving specificity and sensitivity of the clinical tests. In this report, data analyzing will be covered with Python, identifying regions with high and low mutational value as well as identifying single nucleotide polymorphisms (SNPs) in a p53 dataset.
 
.. figure::  https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Figure_workflow.jpg
   :width: 700px
   :height: 450px
   :align: center

   Figure 1: **Workflow for Cyclomics**. A) for library preparation, patient blood is isolated, which contains circulating tumor DNA (ctDNA) fragments. Using a phosphorylated backbone, these fragments can be circularized and subsequently amplified with rolling circle amplification. This leaves very long reads, which are ideal for nanopore sequencing. B) Nanopore sequencing is used to identify nucleotides in a read and after base and variant calling gives fastq and vcf files. After applying a Data filter, true driver mutation can be identified.

.. _1: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _2: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _3: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _4: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _5: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _6: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html
.. _7: http://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/References.html

|
|

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Next.png
   :align:  center
   :width: 100px
   :height: 100px
   :target: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/General%20concept.html