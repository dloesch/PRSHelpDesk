# PRSHelpDesk
 R package for data prep, PRS estimation, and evaluation

## Description 
PRSHelpDesk is a package that assists with data preparation for PRS analyses and contains several functions for evaluating PRS, performing linear combination of multiple PRS, and transforming a PRS to account for genetic ancestry. Also contains a function for estimation of a p+t PRS using a simple window-based pruning approach and a p-value threshold. Finally, includes functions for annotating a PRS with local ancestry information and estimating partial local ancestry specific PRS. Overall, this package was designed for use in cohorts with complex population structure (i.e. admixture).   

## Usage
For data preparation, GWAS summary statistics and target genotype data is loaded using parseGWASSTATS and parseVCF functions. The resulting gData class objects are merged using the mgData constructor function, which performs QC. For PRS estimation, calcPRS() and calcPRS2() constructs a PRS internally or using PLINK, respectively. The weightPRS() function learns weights for the linear combination of PRS. For PRS evaluation, a number of functions are available, including evalPRS() and testPRS(), which differ in their depth. For partial PRS calculated using local ancestry information, annotLAI() annotates a VCF with the LAI and LAiPRS() estimates the partial PRS.  

## Installation
```r
# install.packages("remotes")
remotes::install_github("dloesch/PRSHelpDesk")
```
## Input formats
This package supports BCF, VCF, and PLINK (BED/BIM/FAM) files. PLINK and bcftools need to be installed. 

## References
1. Márquez-Luna C, Loh PR, Price AL. Multi-ethnic polygenic risk scores improve risk prediction in diverse populations. Genet Epidemiol. 2017;41(8):811-823. doi:10.1002/gepi.22083
2. Pain O, Glanville KP, Hagenaars SP, et al. Evaluation of polygenic prediction methodology within a reference-standardized framework. PLOS Genetics. 2021;17(5):e1009021. doi:10.1371/journal.pgen.1009021
3. Marnetto D, Pärna K, Läll K, et al. Ancestry deconvolution and partial polygenic score can improve susceptibility predictions in recently admixed individuals. Nat Commun. 2020;11(1):1628. doi:10.1038/s41467-020-15464-w
4. Lee SH, Goddard ME, Wray NR, Visscher PM. A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012;36(3):214-224.
