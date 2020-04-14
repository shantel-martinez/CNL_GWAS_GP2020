CNL GWAS PHS
================
S.A. Martinez
2020.04.14

**Genome-wide association analysis** of the Cornell Master elite nursery for the Preharvest Sprouting trait using the spike-wetting tests.

> All files in referenced in this document can be downloaded from [CNL\_GWAS\_GP2020](https://github.com/shantel-martinez/CNL_GWAS_GP2020).
> Users will need the following in the same working directory as this CNL\_GWAS\_PHS.Rmd file: [myGD.csv](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/myGD.csv), [myGM.csv](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/myGM.csv), and [allMasterDatawithOHMI\_June2017miss30.RData](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/allMasterDatawithOHMI_June2017miss30.RData)

This analysis follows the [CNL\_Prep\_Stats\_PHS](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/CNL_Prep_Stats_PHS.md) analysis in the [CNL\_GWAS\_GP2020](https://github.com/shantel-martinez/CNL_GWAS_GP2020) project. The environment should contain only the following dataframes:

``` r
load("./CNL_Prep_Stats_PHS.RData")
rm(list= ls()[!(ls() %in% c('PHScGIDblup','PHSwGIDblup','PHSrGIDblup','PHSWhiteComplete','PHSCombComplete','PHSRedComplete'))])
```

If these dataframes are **not** in the current environment, open the [CNL\_Prep\_Stats\_PHS.Rmd](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/CNL_Prep_Stats_PHS.Rmd) file and run all coding chunks. *Be sure to have all the corresponding input files needed to run*
OR Download the [CNL\_Prep\_Stats\_PHS.Rdata](https://github.com/shantel-martinez/CNL_GWAS_GP2020/blob/master/CNL_Prep_Stats_PHS.RData) files from github, save in the working directory this CNL\_GWAS\_PHS.Rmd is in.

PHS Genotype Files
------------------

Genotyping was conducted using genotyping-by-sequencing on the Illumina HiSeq 2000. 11,604 population-based SNP markers were identified after a 30% missing data and a 5% minor allele frequency filter was applied. [N. Santantonio et al., 2019](http://www.genetics.org/content/early/2019/01/24/genetics.118.301851) curated, imputed, and filtered the genotypic data. The `allMasterDatawithOHMI_June2017miss30.RData` file was given to me by N. Santantonio.

Impute the genotpye dataframe `Ximp` and renamed it `myGD`:

``` r
load("../Data Input/allMasterDatawithOHMI_June2017miss30.RData")
# myGD headers: SNP Marker1 Marker2 ... MarkerN
myGD<- Ximp + 1 
rm(pheno, snpInfo, X, Ximp, XnoHet, README)
write.csv(myGD, "./myGD.csv",row.names=TRUE)
```

Ximp is coded as: -1 is homozygous major allele, 0 is heterozygous, and 1 is homozygous minor allele.
`Ximp + 1` changes the code to 0,1,2.

This may not be the most efficent way to add in the taxa names into the first column, but I found that if I write `myGD` as a csv file, then read it back in as `myGD`, it is now a dataframe rather than a large matrix and includes the first column as the GIDs (or genotype names).

``` r
myGD <- read.csv("../Data Input/myGD.csv", head = TRUE,na.string=c(""," ","NA","na","NaN"))
rownames(myGD) = myGD[,1]
myGD[1:4,1:4]
```

    ##                 X S1A_PART1_1197644 S1A_PART1_1238021 S1A_PART1_2991302
    ## cuGS1       cuGS1                 0                 0                 0
    ## cuGS10     cuGS10                 0                 0                 0
    ## cuGS100   cuGS100                 0                 0                 0
    ## cuGS1000 cuGS1000                 0                 0                 0

### Summary of GIDs in study

Summary of GIDs per location with genotypic data

``` r
library("dplyr")
```

    ## Warning: package 'dplyr' was built under R version 3.5.3

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
PHSComb1 <- subset(PHSCombComplete, GID %in% myGD$X)
c<-PHSComb1 %>%
  group_by(Year, Loc) %>%
  summarize(Comb = length(unique(GID)))

PHSWhite1 <- subset(PHSWhiteComplete, GID %in% myGD$X)
w<- PHSWhite1 %>%
  group_by(Year, Loc) %>%
  summarize(White = length(unique(GID)))
PHSred1 <- subset(PHSRedComplete, GID %in% myGD$X)
r<- PHSred1 %>%
  group_by(Year, Loc) %>%
  summarize(Red = length(unique(GID))) 
a <- cbind(c,w,r) 

a %>%
  dplyr::select(Year, Loc, Comb, White, Red)
```

    ## # A tibble: 24 x 5
    ## # Groups:   Year [9]
    ##    Year  Loc      Comb White   Red
    ##    <chr> <chr>   <int> <int> <int>
    ##  1 2010  Helfer    171   116    55
    ##  2 2010  Ketola    208   151    57
    ##  3 2010  Snyder    204   147    57
    ##  4 2011  Helfer    263   190    73
    ##  5 2011  Ketola    266   199    67
    ##  6 2011  McGowen   263   200    63
    ##  7 2012  Helfer    259   180    79
    ##  8 2012  Ketola    266   189    77
    ##  9 2012  Snyder    250   179    71
    ## 10 2013  Helfer     80    51    29
    ## # ... with 14 more rows
