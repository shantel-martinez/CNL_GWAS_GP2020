CNL Data Prep and Statistical Analysis PHS
================
S.A. Martinez
2020.04.14

**Data preparation and statistical analysis** of the Cornell Master elite nursery for the Preharvest Sprouting trait using the spike-wetting tests.

This study used a diversity panel of 1,353 lines from the wheat breeding programs at Cornell University (1,032), Michigan State (112), Ohio State (85), and private companies (51) in addition to landraces (14). The diversity panel also consists of 904 (67%) white kernel color and 449 (33%) were red kernelled. Downstream analyses were conducted on either data subset of white kernel color (white), red kernel color (red), and combined (comb or both) kernel color datasets.

PHS Phenotype Files
-------------------

The panel was planted at three locations per year spanning 2008 to 2016, and two locations in 2017 and 2018. The locations rotated between Helfer, Ketola, McGowen, and Snyder near Ithaca, NY.

For each location, 5 spikes/heads were harvested from each plot at physiological maturity. After five and seven days of dry after-ripening for white and red kernel varieties, respectively, the spikes were misted. **Sprouting scores** were taken after 4 days of misting. The sprouting score was based on a 0-9 scale, where a score of 0 had no visible germination and 9 had 100% visible germination and coleoptile extension (McMaster and Derera, 1987).

``` r
rm(list = ls())
PHSComb <- read.csv("./PHSAll_Comb_20181222.csv", head = TRUE,na.string=c(""," ","NA","na","NaN"), stringsAsFactors=FALSE) 
PHSComb <- PHSComb[!is.na(PHSComb$RawMean),]
```

> > Refer to `PHS Data Exploration.v3.R` to see how I cleaned up and organized the raw phenotype files.

`GID` is a genotype identification number refered to the line tested. The line name is in the `Entry` column.

`Harvest` date at physiological maturity is in Julian Dates, along with the dates `Mist`ed and `Score`d.

`Sprout1`, `Sprout2`,...`Sprout5` refer to the 0-9 score for the 1st spike, 2nd spike, ..., and 5th spike.

The `RawMean` is a raw mean of all 5 sprouting scores.

`Year`, location `Loc`, environment `Env`, days after-ripened `AR`, and `daysMisted` are also included is the data was recorded.

Other agronomic traits may include `PlotYield`, `TW`, `Moisture`, `lodging`, height `Ht`, heading date `HD`, powdery `mildew`, winter hardiness `wh`, and kernel color `KC` if taken.

#### Kernel Color Subsets

Seperating the dataframes into two df: red or white KC based on `$KC` column

``` r
PHSred <- subset(PHSComb, PHSComb$KC=="R" |PHSComb$KC=="r" )
PHSWhite <-subset(PHSComb, PHSComb$KC=="W" |PHSComb$KC=="w" ) 
write.csv(PHSred, "./CNL_Prep_Stats_PHS_files/PHSAll_Red_20191222.csv",row.names=FALSE)
write.csv(PHSWhite, "./CNL_Prep_Stats_PHS_files/PHSAll_White_20191222.csv",row.names=FALSE)

length(unique(PHSWhite$GID)) #White subset
```

    ## [1] 797

``` r
length(unique(PHSred$GID)) #Red subset
```

    ## [1] 333

I can seperate out red vs white `subset(PHSComb, PHSComb$KC=="R" |PHSComb$KC=="r" )`

`rbind(HelRed, KetRed, McGRed, SnyRed)` combined all the red dataframes. Followed by the white and combined dataframes.
Any `GID`s with no kernel color indicated was omitted from downstream analysis in order to prevent mixed kernel color samples affecting the results. *It is still possible there was human error in coding kernel color*

Year needs to be treated as a character rather than an integer.

``` r
PHSred$Year<-as.factor(PHSred$Year)
PHSWhite$Year<-as.factor(PHSWhite$Year)
PHSred$RawMean<-as.numeric(PHSred$RawMean)
PHSWhite$RawMean<-as.numeric(PHSWhite$RawMean)
```

QUALITY CHECK: Red kernels were after-ripened for 7 days while white kernels were after-ripened for 5 days. Omit lines that have conflicting AR time points. However, keep the `NAs` because that is a lack of harvest and misting date recording, not necessarily the wrong AR length.
NOTE: the 8 days AR red kernel color days was one year with a reason. That years environment induced more dormancy than usual, so the sames were after-ripened for one more day longer.

``` r
library(plyr)
count(PHSred, "AR") 
PHSred <- subset(PHSred, PHSred$AR!=5| is.na(PHSred$AR)) 
count(PHSWhite, "AR")
PHSWhite <- subset(PHSWhite, PHSWhite$AR!=7| is.na(PHSWhite$AR)) 
```

The combined 'Comb' dataframe is the red and the white kernel datasets all together.

``` r
PHSComb$Year<-as.character(PHSComb$Year)
PHSComb$RawMean<-as.numeric(PHSComb$RawMean)
length(unique(PHSComb$GID)) #Combined subset
```

    ## [1] 1130

Make GID `###` to GID `cuGS###` by adding a new row, then replace text cuGSOH to just `OH` and cuGSMSU to just `MSU`.

``` r
PHSred$GIDx <- with(PHSred,paste("cuGS",PHSred$GID,sep="")) 
PHSred$GIDx <- gsub("cuGSOH", "OH", PHSred$GIDx)
PHSred$GIDx <- gsub("cuGSMSU", "MSU", PHSred$GIDx)
PHSred$GID <- NULL
names(PHSred)[33] <- "GID"
PHSred$GID[17:25]
```

    ## [1] "cuGS17" "cuGS17" "cuGS18" "cuGS18" "cuGS18" "cuGS19" "cuGS19" "cuGS19"
    ## [9] "cuGS21"

``` r
PHSComb$GIDx <- with(PHSComb,paste("cuGS",PHSComb$GID,sep="")) 
PHSComb$GIDx <- gsub("cuGSOH", "OH", PHSComb$GIDx)
PHSComb$GIDx <- gsub("cuGSMSU", "MSU", PHSComb$GIDx)
PHSComb$GID <- NULL
names(PHSComb)[33] <- "GID"

PHSWhite$GIDx <- with(PHSWhite,paste("cuGS",PHSWhite$GID,sep="")) 
PHSWhite$GIDx <- gsub("cuGSOH", "OH", PHSWhite$GIDx)
PHSWhite$GIDx <- gsub("cuGSMSU", "MSU", PHSWhite$GIDx)
PHSWhite$GID <- NULL
names(PHSWhite)[33] <- "GID"

PHSWhite$GID[4500:4505]
```

    ## [1] "MSU92" "MSU92" "MSU92" "MSU94" "MSU94" "MSU94"

#### Harvest date as factor

``` r
PHSred$Harvest <- as.factor(PHSred$Harvest)  
PHSComb$Harvest <- as.factor(PHSComb$Harvest)  
PHSWhite$Harvest <- as.factor(PHSWhite$Harvest)  
```

#### Harvest x Year

I realized hat it would make sense that just using harvest date as a fixed effect is not the full picture. For example, a harvest date of 161 is not going to have the same effect as a harvest date of 161 the following year, even though they are labeled "161". After bringing this up in lab meeting, we have two routes:

1.  If I still plan to solely using `rrBLUP` for my downstream analysis, I may not be able to keep using `kin.blup` because when you defined the fixed effect terms`fixed = c("col1","col2")`, I havent found a way to define col1xcol2 interaction. Therefore, I need to move to `mixed.solve` or `kinship.blup` for my calculations, and create my own design matrix that includes the harvestxyear interaction.

2.  I can create a new column in the datasets that concatenates the year info and the harvest date info. This column will then be defined as a fixed effect. Note that the column would need to be a character or factor now, since the numeric values have no continuous meaning.

``` r
PHSComb$YrHarv <- paste(PHSComb$Year, PHSComb$Harvest, sep = "")
PHSWhite$YrHarv <- paste(PHSWhite$Year, PHSWhite$Harvest, sep = "")
PHSred$YrHarv <- paste(PHSred$Year, PHSred$Harvest, sep = "")
#Need to remove NA harvest dates
#"2008NA" , "2016NA", "2017NA" , "2011NA" 
PHSComb$YrHarv <- gsub("2008NA" , "NA", PHSComb$YrHarv)
PHSComb$YrHarv <- gsub("2016NA", "NA", PHSComb$YrHarv)
PHSComb$YrHarv <- gsub("2017NA" , "NA", PHSComb$YrHarv)
PHSComb$YrHarv <- gsub("2011NA", "NA", PHSComb$YrHarv)

PHSWhite$YrHarv <- gsub("2008NA" , "NA", PHSWhite$YrHarv)
PHSWhite$YrHarv <- gsub("2016NA", "NA", PHSWhite$YrHarv)
PHSWhite$YrHarv <- gsub("2017NA" , "NA", PHSWhite$YrHarv)
PHSWhite$YrHarv <- gsub("2011NA", "NA", PHSWhite$YrHarv)

PHSred$YrHarv <- gsub("2008NA" , "NA", PHSred$YrHarv)
PHSred$YrHarv <- gsub("2016NA", "NA", PHSred$YrHarv)
PHSred$YrHarv <- gsub("2017NA" , "NA", PHSred$YrHarv)
PHSred$YrHarv <- gsub("2011NA", "NA", PHSred$YrHarv)
```

#### Remove all missing

``` r
library(tidyr)
PHSWhiteComplete <- PHSWhite %>% drop_na(HD,Harvest,YrHarv) #goes from 4558 obs to 3970
PHSRedComplete <- PHSred %>% drop_na(HD,Harvest,YrHarv) #2373 obs to 2226
PHSCombComplete <- PHSComb %>% drop_na(HD,Harvest,YrHarv) #6931 obs to 6196
```

#### Raw Mean Summary

``` r
library(ggplot2)
font<-element_text(face = "bold",  size = 12)
font2<-element_text(size = 12)

ggplot(PHSRedComplete, aes(x = PHSRedComplete$Year, y = PHSRedComplete$RawMean)) +
  geom_boxplot(aes(fill = factor(Loc)), outlier.shape = 1)+ theme_bw() +
  facet_grid(.~Loc)+ theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text = font2,  axis.title = font,strip.text.x = font,plot.title = font)+
  ylab("PHS Raw Mean") + xlab("Environment")+ guides(fill=FALSE) + ylim(0, 9) +
  ggtitle("Red Kernel Color Sprouting")
```

![](CNL_Prep_Stats_PHS_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
#Repeat for PHSWhite and PHSComb
```

![](CNL_Prep_Stats_PHS_files/figure-markdown_github/unnamed-chunk-11-1.png)

Best Fit Model
--------------

### Lowest AIC

Since the data was unbalanced (`PHSWhite`), I can not immeditaely check anova() for models with HD and Harvest, because not all data points have HD or Harvest values, and therefore they are not the same size of dataset.
Which means I needed to subset the data to the values that are complete for HD, Harvest, Env, and so on (`PHSWhiteComplete`). Then reevaluate the models and check for their lower AIC. Is this a valid comparison, if in reality I am using more data for `PHSWhiteVar_3`, versus calculating AIC `anova(PHSWhiteVar_3v2,PHSWhiteVar_4)` with `PHSWhiteVar_3v2` fewer data points?

> Comparing models with [varying samples sizes](https://stats.stackexchange.com/questions/94718/model-comparison-with-aic-based-on-different-sample-size)

After omitting NAs from the dataframes, all models compare:

``` r
library(lme4)
```

    ## Warning: package 'lme4' was built under R version 3.5.3

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
PHSWhiteVar_1 <- lmer(RawMean~(1|GID), data = PHSWhiteComplete) #simple y~x model 
PHSWhiteVar_2 <- lmer(RawMean~(1|GID)+Loc+Year, data = PHSWhiteComplete) 
PHSWhiteVar_3 <- lmer(RawMean~(1|GID)+Env, data = PHSWhiteComplete)
PHSWhiteVar_4 <- lmer(RawMean~(1|GID)+Env+HD+YrHarv, data = PHSWhiteComplete)
```

    ## fixed-effect model matrix is rank deficient so dropping 11 columns / coefficients

``` r
PHSWhiteVar_4b <- lmer(RawMean~(1|GID)+Env+HD, data = PHSWhiteComplete)
PHSWhiteVar_5 <- lmer(RawMean~(1|GID)+Env+HD, data = PHSWhiteComplete)
PHSWhiteVar_7 <- lmer(RawMean~(1|GID)+Env+(1|YrHarv), data = PHSWhiteComplete) #dropping 
PHSWhiteVar_6<- lmer(RawMean~(1|GID)+Env+Year:Harvest+HD, data = PHSWhiteComplete)  #dropping 
```

    ## fixed-effect model matrix is rank deficient so dropping 100 columns / coefficients

``` r
PHSWhiteVar_10 <- lmer(RawMean~(1|GID)+Env+YrHarv, data = PHSWhiteComplete)
```

    ## fixed-effect model matrix is rank deficient so dropping 11 columns / coefficients

``` r
PHSWhiteVar_13 <- lmer(RawMean~(1|GID)+Loc+Year+Year%in%Loc+Year%in%Harvest, data = PHSWhiteComplete)  #dropping 
```

    ## fixed-effect model matrix is rank deficient so dropping 103 columns / coefficients

``` r
PHSWhiteVar_8<- lmer(RawMean~(1|GID)++Loc+Year+Year:Harvest, data = PHSWhiteComplete)  #dropping 
```

    ## fixed-effect model matrix is rank deficient so dropping 88 columns / coefficients

``` r
PHSWhiteVar_11 <- lmer(RawMean~(1|GID)+Loc+(1|Year:Harvest), data = PHSWhiteComplete)  #dropping 
PHSWhiteVar_12 <- lmer(RawMean~(1|GID)+Loc+Year+(1|Year:Harvest), data = PHSWhiteComplete)  #dropping 
PHSWhiteVar_14 <- lmer(RawMean~(1|GID)+Env+(1|Year:Harvest), data = PHSWhiteComplete)  #dropping 
# PHSWhiteVar_15 <- lmer(RawMean~(1|GID)+(1|Year:Harvest), data = PHSWhiteComplete)  #dropping 
PHSWhiteVar_17 <- lmer(RawMean~(1|GID)+(1|Year:Harvest)+Env+HD, data = PHSWhiteComplete)  #dropping

anova(PHSWhiteVar_1,PHSWhiteVar_2)  
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_1: RawMean ~ (1 | GID)
    ## PHSWhiteVar_2: RawMean ~ (1 | GID) + Loc + Year
    ##               Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## PHSWhiteVar_1  3 14681 14699 -7337.3    14675                             
    ## PHSWhiteVar_2 14 13972 14060 -6971.9    13944 730.73     11  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSWhiteVar_3,PHSWhiteVar_4)  
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_3: RawMean ~ (1 | GID) + Env
    ## PHSWhiteVar_4: RawMean ~ (1 | GID) + Env + HD + YrHarv
    ##               Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## PHSWhiteVar_3 26 13808 13971 -6877.8    13756                             
    ## PHSWhiteVar_4 44 13240 13517 -6576.3    13152 603.03     18  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSWhiteVar_7,PHSWhiteVar_5)  
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_7: RawMean ~ (1 | GID) + Env + (1 | YrHarv)
    ## PHSWhiteVar_5: RawMean ~ (1 | GID) + Env + HD
    ##               Df   AIC   BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
    ## PHSWhiteVar_7 27 13369 13539 -6657.6    13315                        
    ## PHSWhiteVar_5 27 13770 13940 -6858.1    13716     0      0          1

``` r
anova(PHSWhiteVar_8,PHSWhiteVar_10) 
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_8: RawMean ~ (1 | GID) + +Loc + Year + Year:Harvest
    ## PHSWhiteVar_10: RawMean ~ (1 | GID) + Env + YrHarv
    ##                Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSWhiteVar_8  34 13342 13555 -6637.0    13274                         
    ## PHSWhiteVar_10 43 13284 13553 -6598.8    13198 76.372      9  8.466e-13
    ##                   
    ## PHSWhiteVar_8     
    ## PHSWhiteVar_10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSWhiteVar_13,PHSWhiteVar_6) 
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_13: RawMean ~ (1 | GID) + Loc + Year + Year %in% Loc + Year %in% 
    ## PHSWhiteVar_13:     Harvest
    ## PHSWhiteVar_6: RawMean ~ (1 | GID) + Env + Year:Harvest + HD
    ##                Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSWhiteVar_13 43 13284 13553 -6598.8    13198                         
    ## PHSWhiteVar_6  44 13240 13517 -6576.3    13152 45.056      1  1.915e-11
    ##                   
    ## PHSWhiteVar_13    
    ## PHSWhiteVar_6  ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSWhiteVar_12,PHSWhiteVar_11)
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_11: RawMean ~ (1 | GID) + Loc + (1 | Year:Harvest)
    ## PHSWhiteVar_12: RawMean ~ (1 | GID) + Loc + Year + (1 | Year:Harvest)
    ##                Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
    ## PHSWhiteVar_11  7 13432 13476 -6708.9    13418                            
    ## PHSWhiteVar_12 15 13423 13517 -6696.3    13393 25.177      8   0.001451 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSWhiteVar_17,PHSWhiteVar_14)
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_14: RawMean ~ (1 | GID) + Env + (1 | Year:Harvest)
    ## PHSWhiteVar_17: RawMean ~ (1 | GID) + (1 | Year:Harvest) + Env + HD
    ##                Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSWhiteVar_14 27 13369 13539 -6657.6    13315                         
    ## PHSWhiteVar_17 28 13326 13502 -6635.1    13270 45.065      1  1.906e-11
    ##                   
    ## PHSWhiteVar_14    
    ## PHSWhiteVar_17 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSWhiteVar_4b,PHSWhiteVar_10) 
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_4b: RawMean ~ (1 | GID) + Env + HD
    ## PHSWhiteVar_10: RawMean ~ (1 | GID) + Env + YrHarv
    ##                Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSWhiteVar_4b 27 13770 13940 -6858.1    13716                         
    ## PHSWhiteVar_10 43 13284 13553 -6598.8    13198 518.58     16  < 2.2e-16
    ##                   
    ## PHSWhiteVar_4b    
    ## PHSWhiteVar_10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSWhiteVar_10,PHSWhiteVar_4)  
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSWhiteComplete
    ## Models:
    ## PHSWhiteVar_10: RawMean ~ (1 | GID) + Env + YrHarv
    ## PHSWhiteVar_4: RawMean ~ (1 | GID) + Env + HD + YrHarv
    ##                Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSWhiteVar_10 43 13284 13553 -6598.8    13198                         
    ## PHSWhiteVar_4  44 13240 13517 -6576.3    13152 45.056      1  1.915e-11
    ##                   
    ## PHSWhiteVar_10    
    ## PHSWhiteVar_4  ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Note: not all combinations are shown, but were tested.

Red Kernel datasets

``` r
PHSRedVar_1 <- lmer(RawMean~(1|GID), data = PHSRedComplete) #simple y~x model 
PHSRedVar_2 <- lmer(RawMean~(1|GID)+Loc+Year, data = PHSRedComplete) 
PHSRedVar_3 <- lmer(RawMean~(1|GID)+Env, data = PHSRedComplete)
PHSRedVar_4 <- lmer(RawMean~(1|GID)+Env+HD+YrHarv, data = PHSRedComplete) 
```

    ## fixed-effect model matrix is rank deficient so dropping 11 columns / coefficients

``` r
PHSRedVar_5 <- lmer(RawMean~(1|GID)+Env+HD, data = PHSRedComplete)
PHSRedVar_10 <- lmer(RawMean~(1|GID)+Env+YrHarv, data = PHSRedComplete) 
```

    ## fixed-effect model matrix is rank deficient so dropping 11 columns / coefficients

``` r
PHSRedVar_6 <- lmer(RawMean~(1|GID)+Env+Year%in%Harvest+HD, data = PHSRedComplete) 
```

    ## fixed-effect model matrix is rank deficient so dropping 101 columns / coefficients

``` r
PHSRedVar_7 <- lmer(RawMean~(1|GID)+Env+HD+Harvest, data = PHSRedComplete) 
PHSRedVar_11 <- lmer(RawMean~(1|GID)+Loc+(1|Year:Harvest), data = PHSRedComplete) 
PHSRedVar_12 <- lmer(RawMean~(1|GID)+Loc+Year+(1|Year:Harvest), data = PHSRedComplete) 
PHSRedVar_13 <- lmer(RawMean~(1|GID)+Loc+Year+Year%in%Loc+Year%in%Harvest, data = PHSRedComplete)
```

    ## fixed-effect model matrix is rank deficient so dropping 104 columns / coefficients

``` r
PHSRedVar_14 <- lmer(RawMean~(1|GID)+Env+(1|Year:Harvest), data = PHSRedComplete) 
PHSRedVar_8 <- lmer(RawMean~(1|GID)+Loc+Year+Year:Harvest, data = PHSRedComplete) 
```

    ## fixed-effect model matrix is rank deficient so dropping 89 columns / coefficients

``` r
PHSRedVar_17 <- lmer(RawMean~(1|GID)+(1|Year:Harvest)+Env+HD, data = PHSRedComplete)  #dropping

anova(PHSRedVar_1,PHSRedVar_2)  #PHSRedVar_2 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_1: RawMean ~ (1 | GID)
    ## PHSRedVar_2: RawMean ~ (1 | GID) + Loc + Year
    ##             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## PHSRedVar_1  3 7799.1 7816.2 -3896.5   7793.1                             
    ## PHSRedVar_2 14 7231.6 7311.5 -3601.8   7203.6 589.46     11  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_3,PHSRedVar_2)  #PHSRedVar_3 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_2: RawMean ~ (1 | GID) + Loc + Year
    ## PHSRedVar_3: RawMean ~ (1 | GID) + Env
    ##             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## PHSRedVar_2 14 7231.6 7311.5 -3601.8   7203.6                             
    ## PHSRedVar_3 26 6995.8 7144.2 -3471.9   6943.8 259.78     12  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_3,PHSRedVar_4)  #PHSRedVar_4 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_3: RawMean ~ (1 | GID) + Env
    ## PHSRedVar_4: RawMean ~ (1 | GID) + Env + HD + YrHarv
    ##             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## PHSRedVar_3 26 6995.8 7144.2 -3471.9   6943.8                             
    ## PHSRedVar_4 43 6659.8 6905.1 -3286.9   6573.8 370.04     17  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_4,PHSRedVar_5)  #PHSRedVar_5 AIC lower 
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_5: RawMean ~ (1 | GID) + Env + HD
    ## PHSRedVar_4: RawMean ~ (1 | GID) + Env + HD + YrHarv
    ##             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## PHSRedVar_5 27 6978.3 7132.4 -3462.2   6924.3                             
    ## PHSRedVar_4 43 6659.8 6905.1 -3286.9   6573.8 350.53     16  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_8,PHSRedVar_10) #PHSRedVar_10 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_8: RawMean ~ (1 | GID) + Loc + Year + Year:Harvest
    ## PHSRedVar_10: RawMean ~ (1 | GID) + Env + YrHarv
    ##              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSRedVar_8  33 6859.2 7047.5 -3396.6   6793.2                         
    ## PHSRedVar_10 42 6694.5 6934.1 -3305.3   6610.5 182.71      9  < 2.2e-16
    ##                 
    ## PHSRedVar_8     
    ## PHSRedVar_10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_7,PHSRedVar_6) #PHSRedVar_10 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_7: RawMean ~ (1 | GID) + Env + HD + Harvest
    ## PHSRedVar_6: RawMean ~ (1 | GID) + Env + Year %in% Harvest + HD
    ##             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
    ## PHSRedVar_7 39 6660.7 6883.2 -3291.3   6582.7                           
    ## PHSRedVar_6 43 6659.8 6905.1 -3286.9   6573.8 8.9062      4    0.06349 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_14,PHSRedVar_12) #PHSRedVar_10 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_12: RawMean ~ (1 | GID) + Loc + Year + (1 | Year:Harvest)
    ## PHSRedVar_14: RawMean ~ (1 | GID) + Env + (1 | Year:Harvest)
    ##              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSRedVar_12 15 6924.4 7010.0 -3447.2   6894.4                         
    ## PHSRedVar_14 27 6762.3 6916.3 -3354.1   6708.3 186.12     12  < 2.2e-16
    ##                 
    ## PHSRedVar_12    
    ## PHSRedVar_14 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_14,PHSRedVar_17) #PHSRedVar_10 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_14: RawMean ~ (1 | GID) + Env + (1 | Year:Harvest)
    ## PHSRedVar_17: RawMean ~ (1 | GID) + (1 | Year:Harvest) + Env + HD
    ##              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSRedVar_14 27 6762.3 6916.3 -3354.1   6708.3                         
    ## PHSRedVar_17 28 6728.6 6888.3 -3336.3   6672.6 35.721      1  2.277e-09
    ##                 
    ## PHSRedVar_14    
    ## PHSRedVar_17 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(PHSRedVar_13,PHSRedVar_11) #PHSRedVar_10 AIC lower
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: PHSRedComplete
    ## Models:
    ## PHSRedVar_11: RawMean ~ (1 | GID) + Loc + (1 | Year:Harvest)
    ## PHSRedVar_13: RawMean ~ (1 | GID) + Loc + Year + Year %in% Loc + Year %in% 
    ## PHSRedVar_13:     Harvest
    ##              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## PHSRedVar_11  7 6944.9 6984.8 -3465.4   6930.9                         
    ## PHSRedVar_13 42 6694.5 6934.1 -3305.3   6610.5 320.36     35  < 2.2e-16
    ##                 
    ## PHSRedVar_11    
    ## PHSRedVar_13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Combined white and red kernel datasets

``` r
PHSCombVar_1 <- lmer(RawMean~(1|GID), data = PHSCombComplete) #simple y~x model 
PHSCombVar_2 <- lmer(RawMean~(1|GID)+Loc+Year, data = PHSCombComplete) 
PHSCombVar_3 <- lmer(RawMean~(1|GID)+Env, data = PHSCombComplete)
PHSCombVar_4 <- lmer(RawMean~(1|GID)+Env+HD+YrHarv, data = PHSCombComplete) 
PHSCombVar_5 <- lmer(RawMean~(1|GID)+Env+HD, data = PHSCombComplete)
PHSCombVar_10 <- lmer(RawMean~(1|GID)+Env+YrHarv, data = PHSCombComplete) 
PHSCombVar_6 <- lmer(RawMean~(1|GID)+Env+(1|Year:Harvest), data = PHSCombComplete) 
PHSCombVar_7 <- lmer(RawMean~(1|GID)+Env+HD+Harvest, data = PHSCombComplete) 
PHSCombVar_8 <- lmer(RawMean~(1|GID)+Loc+Year+Year:Harvest, data = PHSCombComplete) 

PHSCombVar_15 <- lmer(RawMean~(1|GID)+Env+KC, data = PHSCombComplete)

PHSCombVar_14 <- lmer(RawMean~(1|GID)+Env+(1|Year:Harvest), data = PHSCombComplete)
PHSCombVar_11 <- lmer(RawMean~(1|GID)+Loc+(1|Year:Harvest), data = PHSCombComplete) 
PHSCombVar_12 <- lmer(RawMean~(1|GID)+Loc+Year+(1|Year:Harvest), data = PHSCombComplete) 
PHSCombVar_13 <- lmer(RawMean~(1|GID)+Loc+Year+Year%in%Loc+(1|Year:Harvest), data = PHSCombComplete) 
PHSCombVar_16 <- lmer(RawMean~(1|GID)+Env+(1|Year:Harvest)+HD+KC, data = PHSCombComplete)
PHSCombVar_17 <- lmer(RawMean~(1|GID)+Env+(1|Year:Harvest)+HD, data = PHSCombComplete)

# anova(PHSCombVar_1,PHSCombVar_2)  #PHSCombVar_2 AIC lower
# anova(PHSCombVar_3,PHSCombVar_2)  #PHSCombVar_3 AIC lower
# anova(PHSCombVar_3,PHSCombVar_4)  #PHSCombVar_3 AIC lower
# anova(PHSCombVar_15,PHSCombVar_8)  #PHSCombVar_8 AIC lower
# anova(PHSCombVar_3,PHSCombVar_5)  #PHSCombVar_5 AIC lower
# anova(PHSCombVar_5,PHSCombVar_10) #PHSCombVar_10 AIC lower
# anova(PHSCombVar_7,PHSCombVar_6)  #PHSCombVar_10 AIC lower
# anova(PHSCombVar_17,PHSCombVar_11) #PHSCombVar_10 AIC lower
# anova(PHSCombVar_13,PHSCombVar_12) #PHSCombVar_10 AIC lower
# anova(PHSCombVar_14,PHSCombVar_16) #PHSCombVar_10 AIC lower
```

For all of the models with a Year x harvest interation stated in the equation, there ae many values/columns dropped due to insufficient data. This could be because there was only one or two data points per conditions (i.e., `181_2008` n = 2). By using a concatenated fixed effect (no interaction defined), fewer columns are dropped. The difference of the AIC values between the interaction term `Year:Harvest` and the concatenated `YrHarv` do not appear to be drastic.

Check heritability to see if using either effect in the model changes h2.

### Calculating Heritability

Basic heritability is h2 = Vg / (Vg + Ve), where Vg is the genetic variance `GID Variance` and Ve is the `residual variance` (i.e. variance due to non genetics effects such as environment)

However, since this is an unbalanced dataset, and I'm eventually going to calculate BLUPs, the generalized measure of heritability should be used from [Cullis et al., 2006](https://doi.org/10.1198/108571106X154443)

``` r
Cullis_H2=function(model){
  library(arm)
  ses<- se.ranef(model)$'GID'
  v_BLUP<- ses^2
  sigma2_g=VarCorr(model, comp="Variance")$'GID'[1]
  Reliability<- 1- v_BLUP/ (2*sigma2_g) 
  H2<- round(mean(Reliability),3)
  H2
}
```

where `model` is your defined model object from `lmer`
where `sigma2_g` is the genetic variance estimated with the model saved in `model`
`H2` is equivalent to broad-sense heritability on the line-mean (or family-mean, if your individuals are non-inbred families) basis

Combined dataset h2

``` r
h2_1 <- Cullis_H2(PHSCombVar_1)
h2_2 <- Cullis_H2(PHSCombVar_2)
h2_3 <- Cullis_H2(PHSCombVar_3)
h2_15 <- Cullis_H2(PHSCombVar_15)
h2_4 <- Cullis_H2(PHSCombVar_4)
h2_10 <- Cullis_H2(PHSCombVar_10)
h2_8 <- Cullis_H2(PHSCombVar_8)
h2_14 <- Cullis_H2(PHSCombVar_14)
h2_16 <- Cullis_H2(PHSCombVar_16)
h2_12 <- Cullis_H2(PHSCombVar_12)
h2_11 <- Cullis_H2(PHSCombVar_11)
h2_17 <- Cullis_H2(PHSCombVar_17)

h2 <- c(h2_1, h2_2, h2_3,h2_15, h2_4, h2_10, h2_8,h2_14,h2_16, h2_12,h2_11,h2_17)
print("Comb h2");sprintf('%.2f',h2)
```

    ## [1] "Comb h2"

    ##  [1] "0.87" "0.88" "0.88" "0.80" "0.89" "0.89" "0.89" "0.89" "0.82" "0.89"
    ## [11] "0.89" "0.89"

`RawMean~(1|GID)+Env+HD+YrHarv` gets h2 = 0.89; note that if I use this mixed model, I will lose some data that has HD and Harvest missing. Every time I add KC into my mixed model, h2 drastically reduces (0.89 to 0.80). The residual is much higher compared to the genetic variance. This is even after combing through the data for unclear kernel color observations and omitting any conflicting data points.

Repeat for white and red kernel color datasets

    ## [1] "White h2 : h2_1, h2_2, h2_3, h2_4, h2_10, h2_6, h2_8,h2_14, h2_12,h2_11,h2_17"

    ##  [1] "0.77" "0.80" "0.81" "0.82" "0.83" "0.82" "0.83" "0.83" "0.83" "0.83"
    ## [11] "0.82"

    ## [1] "Red h2 : h2_1, h2_2, h2_3, h2_4, h2_10, h2_6, h2_8,h2_14, h2_12,h2_11,h2_17"

    ##  [1] "0.80" "0.83" "0.85" "0.86" "0.87" "0.86" "0.86" "0.87" "0.86" "0.86"
    ## [11] "0.86"

`RawMean~(1|GID)+Env+HD+YrHarv` gets h2 = 0.82 for white KC and 0.86 for red KC; note that if I use this mixed model, I will lose a lot of data that has HD and Harvest missing, however the model fit is better with increased h2.

### BLUP Analysis

Once we have the model covariates we wish to use `RawMean~(1|GID)+Env+HD+YrHarv`, then we can calculate the BLUPs using that model.

``` r
blup <- ranef(PHSWhiteVar_4,drop=TRUE)
y <- attributes(blup[[1]])

blup_W <- data.frame(blup[["GID"]])
rownames(blup_W) <- y$`names`
blup_W$GID <- rownames(blup_W)


blup <- ranef(PHSRedVar_4,drop=TRUE)
y <- attributes(blup[[1]])

blup_r <- data.frame(blup[["GID"]])
rownames(blup_r) <- y$names
blup_r$GID <- rownames(blup_r)

blup <- ranef(PHSCombVar_4,drop=TRUE)
y <- attributes(blup[[1]])

blup_C <- data.frame(blup[["GID"]])
rownames(blup_C) <- y$names
blup_C$GID <- rownames(blup_C)
```

Now `PHScGIDblup` `PHSwGIDblup` and `PHSrGIDblup` are the dataframes with BLUP values in col \[1\] and the GID names in the row names.
However, we will be comparing Comb, Red, and White df to one another, therefore a new column `KC` needs to be defined to differentiate the three.

``` r
PHSwGIDblup <-  blup_W
PHScGIDblup <-  blup_C
PHSrGIDblup <-  blup_r
PHScGIDblup$KC <- "Comb"
PHSrGIDblup$KC <- "Red"
PHSwGIDblup$KC <- "White"
PHScGIDblup$GID <- row.names(PHScGIDblup)
PHSrGIDblup$GID <- row.names(PHSrGIDblup)
PHSwGIDblup$GID <- row.names(PHSwGIDblup)
```

### BLUP and Raw plot

Plot the raw means within each environment and the BLUPs across all environments.

``` r
a <- ggplot(PHSWhiteComplete, aes(x = PHSWhiteComplete$Year, y = PHSWhiteComplete$RawMean))+
  geom_boxplot(aes(fill = factor(Loc)), outlier.shape = 1)+ theme_bw() +
  facet_grid(.~Loc)+ theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.4))+
  theme(axis.text = font2,  axis.title = font,strip.text.x = font,plot.title = font)+
  ylab("Sprout Raw Mean") + xlab("Environment")+ guides(fill=FALSE) + ylim(0, 9) +
  ggtitle("A)")+scale_fill_grey()

PHSwGIDblup$Env <- "All Env"
b <- ggplot(PHSwGIDblup, aes(PHSwGIDblup$KC, PHSwGIDblup$blup...GID...)) +
  geom_boxplot(aes(fill = factor(KC)))+ theme_bw() +
    scale_fill_manual(values="#FFD39B")+
  theme(axis.text = font2,  axis.title = font, strip.text.x = font, 
        plot.title = font, legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.4))+
  ylab("Sprout BLUP") + xlab("") + ylim(-3,3)+
  facet_grid(.~Env)+ggtitle("B)")+
  annotate("text", x =1, y = -2.5, label = paste("italic(h) ^ 2 == ",h2w),  parse = TRUE,size =6)

c <- ggplot(PHSRedComplete, aes(x = PHSRedComplete$Year, y = PHSRedComplete$RawMean))+
  geom_boxplot(aes(fill = factor(Loc)), outlier.shape = 1)+ theme_bw() +
  facet_grid(.~Loc)+ theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.4))+
  theme(axis.text = font2,  axis.title = font,strip.text.x = font,plot.title = font)+
  ylab("Sprout Raw Mean") + xlab("Environment")+ guides(fill=FALSE) + ylim(0, 9) +
  ggtitle("C)")+scale_fill_grey()

PHSrGIDblup$Env <- "All Env"
d <- ggplot(PHSrGIDblup, aes(PHSrGIDblup$KC, PHSrGIDblup$blup...GID...)) +
  geom_boxplot(aes(fill = factor(KC)))+ theme_bw() +
    scale_fill_manual(values="#8B0000")+
  theme(axis.text = font2,  axis.title = font, strip.text.x = font, 
        plot.title = font, legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.4))+
  ylab("Sprout BLUP") + xlab("") + ylim(-3,3)+
  facet_grid(.~Env)+ggtitle("D)")+
  annotate("text", x =1, y = -2.5, label = paste("italic(h) ^ 2 == ",h2r),  parse = TRUE,size =6)

library(gridExtra)
a <- a+ xlab("")
p <- grid.arrange(a, b, c,d, nrow = 2,  widths = c(5, 1),heights = c(14,12))
```

![](CNL_Prep_Stats_PHS_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
# ggsave("./CNL_Prep_Stats_PHS_files/SproutMeans20200409.png", plot =p,width = 16, height = 12, units = "in")
```

Before we move any further to downstream analysis, we've accumulated a lot of dataframes in our environment while prepping the three BLUP dataframes.

Tidy up the global environment to only include `PHScGIDblup`, `PHSwGIDblup` , and `PHSrGIDblup` .

``` r
rm(list= ls()[!(ls() %in% c('PHScGIDblup','PHSwGIDblup','PHSrGIDblup','PHSWhiteComplete','PHSCombComplete','PHSRedComplete'))])
```
