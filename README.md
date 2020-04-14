# CNL_GWAS_GP2020  
Authors: SA Martinez, DW Sweeney, ME Sorrells
v.2020.04.14

Complete analysis of the Cornell Master elite nursery for the Preharvest Sprouting trait using the spike-wetting tests.  

This study used a diversity panel of 1,353 lines from the wheat breeding programs at Cornell University (1,032), Michigan State (112), Ohio State (85), and private companies (51) in addition to landraces (14). The diversity panel also consists of 904 (67%) white kernel color and 449 (33%) were red kernelled. Downstream analyses were conducted on either data subset of white kernel color (white) or red kernel color (red) datasets.  

## Summary of Analysis Files

| Description          | .Rmd File Name          | Input Files Needed         |
| --------------------- | ---------------- | ------------------------- |
| Data Prep, BLUP, h2    | CNL_Prep_Stats_PHS.Rmd | PHSAll_Comb_20181222.csv |
| Genome-wide association analysis      | CNL_GWAS_PHS.Rmd | myGD.csv |
|        |   | myGM.csv |
|        |   | allMasterDatawithOHMI_June2017miss30.RData |
| Genomic prediction - all environments |            |         |

All input files to reproduce the analysis can be downloaded from [CNL_GWAS_GP2020](https://github.com/shantel-martinez/CNL_GWAS_GP2020/tree/master/Data%20Input).      

## Data Prep and Statistics  
Viewable data analysis: [CNL_Prep_Stats_PHS](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/Data%20Analysis/CNL_Prep_Stats_PHS.md)   
Downloadable Rmd File: [CNL_Prep_Stats_PHS.Rmd](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/Data%20Analysis/CNL_Prep_Stats_PHS.Rmd) 

## GWAS  
Viewable data analysis: [CNL_GWAS_PHS](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/Data%20Analysis/CNL_GWAS_PHS.md) 
Downloadable Rmd File: [CNL_GWAS_PHS.Rmd](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/Data%20Analysis/CNL_GWAS_PHS.Rmd)  

