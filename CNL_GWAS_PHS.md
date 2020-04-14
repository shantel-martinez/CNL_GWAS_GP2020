CNL GWAS PHS
================
S.A. Martinez
2020.04.14

**Genome-wid association analysis** of the Cornell Master elite nursery for the Preharvest Sprouting trait using the spike-wetting tests.

This study used a diversity panel of 1,353 lines from the wheat breeding programs at Cornell University (1,032), Michigan State (112), Ohio State (85), and private companies (51) in addition to landraces (14). The diversity panel also consists of 904 (67%) white kernel color and 449 (33%) were red kernelled. Downstream analyses were conducted on either data subset of white kernel color (white), red kernel color (red), and combined (comb or both) kernel color datasets.

> All files in referenced in this document can be downloaded from [CNL\_GWAS\_GP2020](https://github.com/shantel-martinez/CNL_GWAS_GP2020).
> Users will need the following in the same working directory as this CNL\_Prep\_Stats\_PHS.Rmd file: [myGD.csv](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/myGD.csv), [myGM.csv](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/myGM.csv), and [allMasterDatawithOHMI\_June2017miss30.RData](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/allMasterDatawithOHMI_June2017miss30.RData)

This analysis follows the [CNL\_Prep\_Stats\_PHS](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/CNL_Prep_Stats_PHS.md) analysis in the [CNL\_GWAS\_GP2020](https://github.com/shantel-martinez/CNL_GWAS_GP2020) project. The environment should contain only the following:

``` r
rm(list= ls()[!(ls() %in% c('PHScGIDblup','PHSwGIDblup','PHSrGIDblup','PHSWhiteComplete','PHSCombComplete','PHSRedComplete'))])
```

If these dataframes are **not** in the current environment, open the [CNL\_Prep\_Stats\_PHS.Rmd](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/CNL_Prep_Stats_PHS.Rmd) file and run all coding chunks.
