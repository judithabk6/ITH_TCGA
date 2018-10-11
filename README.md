# ITH_TCGA

this repository contains the code necessary to run the analysis of ITH (intra-tumor heterogeneity) in 3 cancer types on TCGA data published in "Assessing reliability of intra-tumor heterogeneity estimates from single sample whole exome sequencing data" Abecassis et al, 2018. 
The master code can be found in the bash script ```run_all.sh```. It is not advised to execute it all at once, considering the run time, and the paths to change (python virtual environments, paths to executables etc).

#### Installation and dependencies
this project uses many packages. Here is a list of the most important ones. Requirement files are provided to use with pip. Here is a list of the main ones, including some that can be installed with pip (specified). A Python or R version is specified when it did not work with R-3.2.3 or Python 3.5.5. For Python, 3 virtual environments were used for the 3 needed versions (2.7, 3.5, 3.6). All computation was performed under a Debian distribution, with torque as scheduler.

##### data download packages
- [gdc transfer tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
- [TCGABiolinks 2.9.5](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) 	(R >= 3.3)

##### data preprocessing
- [segment_liftover 0.951](https://pypi.org/project/segment-liftover/) (python >= 3.6 - installed by pip)
- [DESeq2 2_1.14.0](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

##### ITH methods
- [PyClone 0.13.0](https://bitbucket.org/aroth85/pyclone/wiki/Installation) (Python 2.7)
- [PhyloWGS 8b99d71ec0](https://github.com/morrislab/phylowgs/tree/8b99d71ec05d55dfb6f5bf914eac4817cb2731d5) (Python 2.7)
- [Sciclone 1.1](https://github.com/genome/sciclone)
- [EXPANDS 2.0.0](https://cran.r-project.org/web/packages/expands/index.html)
- [CSR 901f4c7](https://github.com/kaixiany/CSR/tree/901f4c7418373e6de5835607f3fbb8c53bc7128e)

##### statistical packages
- [survcomp 1.26.0](https://bioconductor.org/packages/release/bioc/html/survcomp.html) (R >= 3.4)
- [scikit-survival](https://github.com/sebp/scikit-survival)

All scripts should be run from the root of the folder. Folder tmp contains the final file with all results from ITH methods, allowing to reproduce the survival analysis.
