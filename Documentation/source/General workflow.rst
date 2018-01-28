General workflow
----------------
**Python programming**. First, python language understanding and writing was conducted by using several test exercises not incorporated into this report. All python codes were made into jupyter notebook from anaconda. Python programming and visualization included several library packages, such as pandas, numpy and bokeh. All functions and scripts in Nanopore Data Analysis (NaDA_) were self written and involve only python as language. 

**Data analysis**. Data analysis was divided between fastq and vcf files. Goals for these files were respectively finding high quality score sequences and identifying mutation present in a p53 wild-type database in comparison to a p53 mutated database. All analysis on these files was done using python by self written scripts. Overarching goal for both these files were to acquire knowledge on ctDNA sequencing in order to create a data filter. This report doesn't cover actual filtering as data analysis could not conclude any specific regions.

**High-Performance Computing (HPC)**. Running scripts on larger datasets require much more memory and an overall better computer performance. Therefore, the HPC facility was used for data analysis. The HPC facility consists of 1544 cores and 600TB storage for fast and parallel data analysis. Currently, project NaDA_ is easily accessible on the HPC for further analysis and visualization of Nanopore datasets.

.. _NaDA: https://github.com/DouweSpaanderman/NaDA/

|
|

.. figure:: https://raw.githubusercontent.com/DouweSpaanderman/NaDA/master/Documentation/source/_static/Next.png
   :align: center
   :target: https://htmlpreview.github.io/?https://github.com/DouweSpaanderman/NaDA/blob/master/Documentation/build/html/Results%20and%20Discussion.html