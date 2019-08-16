---
title: "Labcyte-RT Data Analysis (8th experiment)"
subtitle: "RT optimisation with the Labcyte Echo 525: TSO, RT primer and RNA amounts"
output:
  html_document:
    keep_md: yes
    number_sections: yes
    toc: yes
    pandoc_args:
     - "--lua-filter=links-to-markdown.lua"
editor_options: 
  chunk_output_type: console
---



Experiment 8 follows exactly the same design as for
[experiment 7](Labcyte-RT_Data_Analysis_7.md), except that we are using
SuperScript IV and its buffer.

In this file, we analyse the results for SSIII and SSIV in parallel.

Here, we assessed again multiple combinations of TSO, RT primer and RNA amounts,
using a different stock of TSOs (PO_8268526), purchased earlier but apparently
of better quality (see [experiment 6](Labcyte-RT_Data_Analysis_6.md)), and
with a more extensive randomisation of TSO barcodes and well coordinates
(see designs [6a](Labcyte-RT6a.md), [6b](Labcyte-RT6b.md), [6c](Labcyte-RT6c.md) and [6d](Labcyte-RT6d.md)).


Load R packages
===============


```r
library("CAGEr")
library("ggplot2")
library("magrittr")
library("MultiAssayExperiment")
library("SummarizedExperiment")
library("viridis")
```


Load CAGE libraries
===================


```r
if(file.exists("Labcyte-RT_Data_Analysis_7+8.Rds")) {
  ce <- readRDS("Labcyte-RT_Data_Analysis_7+8.Rds")
} else {
  ce7 <- readRDS(paste0("Labcyte-RT_Data_Analysis_", 7, ".Rds"))
  ce8 <- readRDS(paste0("Labcyte-RT_Data_Analysis_", expNumber, ".Rds"))
  ce <- mergeCAGEsets(ce7, ce8)
  ce$enzyme <- ce$plateID
  levels(ce$enzyme) <- c(rep("SSIII", 4), rep("SSIV", 4))
  ce$group <- paste0(ce$RNA_factor, ", ", ce$enzyme) %>% factor(
    c( "100 ng RNA, SSIII", "100 ng RNA, SSIV"
     ,  "10 ng RNA, SSIII",  "10 ng RNA, SSIV"
     ,   "1 ng RNA, SSIII",   "1 ng RNA, SSIV"
     , "100 pg RNA, SSIII", "100 pg RNA, SSIV"
     ,  "10 pg RNA, SSIII",  "10 pg RNA, SSIV"
     ,   "1 pg RNA, SSIII",   "1 pg RNA, SSIV"
     ))
  ce$group.lattice <- paste0(ce$RNA_factor, ", ", ce$enzyme) %>% factor(
    c(   "1 pg RNA, SSIII",   "1 pg RNA, SSIV"
     ,  "10 pg RNA, SSIII",  "10 pg RNA, SSIV"
     , "100 pg RNA, SSIII", "100 pg RNA, SSIV"
     ,   "1 ng RNA, SSIII",   "1 ng RNA, SSIV"
     ,  "10 ng RNA, SSIII",  "10 ng RNA, SSIV"
     , "100 ng RNA, SSIII", "100 ng RNA, SSIV"
     )) # ggplot2 and lattice use a different order...
  ce$group_horiz <- paste0(ce$RNA_factor, ", ", ce$enzyme) %>% factor(
    c(   "1 pg RNA, SSIII",  "10 pg RNA, SSIII" , "100 pg RNA, SSIII"
     ,   "1 ng RNA, SSIII",  "10 ng RNA, SSIII" , "100 ng RNA, SSIII"
     ,   "1 pg RNA, SSIV" ,  "10 pg RNA, SSIV"  , "100 pg RNA, SSIV"
     ,   "1 ng RNA, SSIV" ,  "10 ng RNA, SSIV"  , "100 ng RNA, SSIV"
     ))
  annotateCTSS(ce, rtracklayer::import.gff("gencode.vM1.annotation.gtf.gz"))
  saveRDS(ce, "Labcyte-RT_Data_Analysis_7+8.Rds")
}
```


Remove negative controls
========================


```r
ce <- ce[, ce$RNA_vol    != 0]
ce <- ce[, ce$RT_PRIMERS != 0]
```


Analysis
========

## Color code

In legends and axis panels:

 - RNA amounts are written in blue;
 - RT primer molarities are in green;
 - TSO molarities are in brown/red.
 
Discrepancies signal a mislabeling.


```r
theme_TSO_by_RTP_facet_RNA <- function() {
  theme( axis.text.x      = element_text(colour = "brown")
       , axis.title.x     = element_text(colour = "brown")
       , legend.title     = element_text(colour = "darkgreen")
       , legend.text      = element_text(colour = "darkgreen")
       , strip.background = element_rect(fill   = NA)
       , strip.text       = element_text(colour = "darkblue"))
}

theme_TSO_by_RNA_facet_RNA <- function() {
  theme( axis.text.x      = element_text(colour = "brown")
       , axis.title.x     = element_text(colour = "brown")
       , legend.title     = element_text(colour = "darkblue")
       , legend.text      = element_text(colour = "darkblue")
       , strip.background = element_rect(fill   = NA)
       , strip.text       = element_text(colour = "darkblue"))
}

theme_RTP_by_TSO_facet_RNA <- function() {
  theme( axis.text.x      = element_text(colour = "darkgreen")
       , axis.title.x     = element_text(colour = "darkgreen")
       , legend.title     = element_text(colour = "brown")
       , legend.text      = element_text(colour = "brown")
       , strip.background = element_rect(fill   = NA)
       , strip.text       = element_text(colour = "darkblue"))
}

theme_RTP_by_RNA_facet_RNA <- function() {
  theme( axis.text.x      = element_text(colour = "darkgreen")
       , axis.title.x     = element_text(colour = "darkgreen")
       , legend.title     = element_text(colour = "darkblue")
       , legend.text      = element_text(colour = "darkblue")
       , strip.background = element_rect(fill   = NA)
       , strip.text       = element_text(colour = "darkblue"))
}
```

### QC indicators


```r
ce$TD_rate <- ce$tagdust / ce$extracted * 100
ce$rRNA_rate <- ce$rdna / (ce$extracted - ce$tagdust) * 100
ce$mapping_rate <- ce$properpairs / (ce$extracted - ce$tagdust) * 100
ce$strand_invasion_rate <- ce$strandInvaders / (ce$counts + ce$strandInvaders) * 100
ce$promoter_rate <- ce$promoter / (ce$counts) * 100
```

#### Richness


```r
CTSStoGenes(ce)
ce$r10g <- vegan::rarefy(t(assay(GeneExpSE(ce))), 10)
```

```
## Warning in vegan::rarefy(t(assay(GeneExpSE(ce))), 10): requested 'sample'
## was larger than smallest site maximum (0)
```

```r
ce$r10g[ce$counts < 10] <- NA
ce$r100g <- vegan::rarefy(t(assay(ce[["geneExpMatrix"]])),100)
```

```
## Warning in vegan::rarefy(t(assay(ce[["geneExpMatrix"]])), 100): requested
## 'sample' was larger than smallest site maximum (0)
```

```r
ce$r100g[ce$counts < 100] <- NA 
```

#### Yield

Because we multiplexed reactions together, the ones with the highest yield
will give the largest amount of reads.  Higher yield gives the possibility
of reducing the number of PCR cycles.

Since multiplexing is not perfect, each library had a different number
of reads.  Therefore, to compare yields in terms of number of aligned reads,
etc, one needs to normalise per indexed library.


```r
tapply(ce$librarySizes, paste(ce$index, ce$plateID), sum)
```

```
##  AAGAGGCA R AAGAGGCA R2  ACTCGCTA S ACTCGCTA S2  ACTGAGCG T ACTGAGCG T2 
##      123337       53378      127701       73249      131305       54776 
##  AGGCAGAA Q AGGCAGAA Q2  ATCTCAGG R ATCTCAGG R2  ATGCGCAG S ATGCGCAG S2 
##      115271       54025       20575       11084        3363       10874 
##  CCTAAGAC T CCTAAGAC T2  CGAGGCTG R CGAGGCTG R2  CGATCAGT T CGATCAGT T2 
##       94921       59697      141531       77501       76723       54713 
##  CGGAGCCT S CGGAGCCT S2  CGTACTAG Q CGTACTAG Q2  CTCTCTAC R CTCTCTAC R2 
##      123633       53897      104142       72227      228370       92665 
##  GCGTAGTA S GCGTAGTA S2  GCTCATGA R GCTCATGA R2  GGACTCCT Q GGACTCCT Q2 
##      121044       80485       90912       23407       33080       35209 
##  GGAGCTAC S GGAGCTAC S2  GTAGAGGA R GTAGAGGA R2  TAAGGCGA Q TAAGGCGA Q2 
##      116273       70829      134703       56759      129884       83945 
##  TACGCTGC S TACGCTGC S2  TAGCGCTC T TAGCGCTC T2  TAGGCATG Q TAGGCATG Q2 
##       27205       28943      127514      119531        2063        9414 
##  TCCTGAGC Q TCCTGAGC Q2  TCGACGTC T TCGACGTC T2  TGCAGCTA T TGCAGCTA T2 
##       63777       40199        3490       10377       27828       21581
```

```r
libMean <- tapply(ce$librarySizes, paste(ce$index, ce$plateID), mean)

ce$libSizeNormBylib <-
  mapply( FUN   = function(n, index, plateID)
                    n / libMean[paste(index, plateID)]
        , n     = ce$librarySizes
        , index = ce$index
        , plateID = ce$plateID)
```

#### Aggregated table for contour plots


```r
z <- aggregate(colData(ce)[,c("TD_rate", "rRNA_rate", "libSizeNormBylib", "mapping_rate", "strand_invasion_rate", "promoter_rate", "r10g", "r100g")], by=list(TSO=ce$TSO, RTP=ce$RT_PRIMERS, group=ce$group), median) %>% as.data.frame()
z$group.lattice <- factor(z$group, levels=levels(ce$group.lattice))

zm <- aggregate(colData(ce)[,c("TD_rate", "rRNA_rate", "libSizeNormBylib", "mapping_rate", "strand_invasion_rate", "promoter_rate", "r10g", "r100g")], by=list(TSO=ce$TSO, RTP=ce$RT_PRIMERS, group=ce$group), function(x) mad(x) / median(x)) %>% as.data.frame()
zm$group.lattice <- factor(z$group, levels=levels(ce$group.lattice))

zm$avg <- rowMeans(cbind(zm$TD_rate, zm$rRNA_rate, zm$mapping_rate, zm$strand_invasion_rate, zm$promoter_rate))

# Calculated mads to survey variability but the plots are uninspiring.

asp <- diff(range(log(z$RTP))) / diff(range(log(z$TSO))) # for lattice::contourplot
```

#### Lattice helper functions


```r
lowerIsBetter <- function(data, zvalue, main, cuts = 10) {
  lattice::contourplot(
     data = data,
     zvalue ~ TSO * RTP | group.lattice,
     cuts = cuts,
     scales=list(x=list(at= c(0.6, 2.5, 10, 40, 160), log=10),
                 y=list(at=c(1, 4, 16), log=10)),
     aspect = asp,
     main = main,
     xlab = "Template-switching oligonucleotides (μM)",
     ylab = "Reverse transcription primers (μM)",
     region = TRUE,
     col.regions = colorRampPalette(c("cyan", "white", "pink"))(100))
}

higherIsBetter <- function(data, zvalue, main, cuts = 10) {
  lattice::contourplot(
     data = data,
     zvalue ~ TSO * RTP | group.lattice,
     cuts = cuts,
     scales=list(x=list(at= c(0.6, 2.5, 10, 40, 160), log=10),
                 y=list(at=c(1, 4, 16), log=10)),
     aspect = asp,
     main = main,
     xlab = "Template-switching oligonucleotides (μM)",
     ylab = "Reverse transcription primers (μM)",
     region = TRUE,
     col.regions = colorRampPalette(c("pink", "white", "cyan"))(100))
}
```

## Tag dust

If not removed during library preparation, oligonucleotide artifacts strongly
dominate libraries prepared with 1 pg RNA.  In general, the amount of artefacts
increases when starting RNA amounts decrease.  Here, increasing RT primer
molaritys increase artefacts.  In contrary, and somewhat surprisingly,
increasing TSO molarity seems to reduce artefacts.

Sub panel at 10,000 pg is noisy because replicate `CGAGGCTG` is an outlier with
a large amount of artefacts.


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, TD_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10.0, 40.0, 160.0)) +
  scale_y_continuous("Tag dust (% of extracted reads)") +
  labs( title = "Amount of oligonucleotide artefacts"
      , caption = stringr::str_wrap(width = 50, "The amount of artefacts detected by TagDust is increased by RT primers and decreased by RNA and TSOs.  SuperScript IV generated less artefacts than SSIII.")) +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_tagDust_by_RTP_facet_RNA-1.png)<!-- -->

Sub panel at 10,000 pg is noisy because replicate `CGAGGCTG` is an outlier with
a large amount of artefacts.  Here is the same figure with that library removed.


```r
ggplot(colData(ce[, ! (ce$index == "CGAGGCTG" & ce$plateID == "R")]) %>% data.frame, aes(TSO, TD_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10.0, 40.0, 160.0)) +
  scale_y_continuous("Tag dust (% of extracted reads)") +
  labs( title = "Amount of oligonucleotide artefacts"
      , subtitle = '(outlier library "CGAGGCTG" of plate "R" removed)'
      , caption = stringr::str_wrap(width = 50, "The amount of artefacts detected by TagDust is increased by RT primers and decreased by RNA and TSOs.  SuperScript IV generated less artefacts than SSIII.")) +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_tagDust_by_RTP_facet_RNA_no_outlier-1.png)<!-- -->

### Contour plot


```r
lowerIsBetter(z, z$TD_rate, "Amount of oligonucleotide artefacts (%)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TagDust_countour-1.png)<!-- -->

## Ribosomal RNA

 - RT primers molarity increases rRNA rate.
 - With SSIII, rRNA rate is maximal at mild amounts of RNA (~1 ng) and at
   high TSO concentration, a minimum is reached between 20 and 40 µM,
   depending on the quantity RNA and RT primers.


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, rRNA_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("rRNA (% of non-tagdust extracted reads)") +
  labs( title = "Fraction of reads aligning to rRNA sequences."
      , caption = stringr::str_wrap(width = 50, "The fraction of reads aligning to ribosomal RNA (rRNA) sequences (after removing oligonucleotide artefacts) is increased by RT primers.  It varies with TSO and RNA amounts.  Overall, SuperScript IV gives more rRNA reads than SSIII.")) +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_rRNA_by_RTP_facet_RNA-1.png)<!-- -->

### Contour plot


```r
lowerIsBetter(z, z$rRNA_rate, "Amount of ribosomal RNA sequences (%)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/rRNA_countour-1.png)<!-- -->


## Yield

RT primer molarity mildly influences yield.  Higher molarities are
needed when TSO molarity is increased.  Conversely, high molarities are
detrimental for low TSO amounts.  In brief, the RT primer concentration must
be adjusted to the TSO concentration.


```r
ggplot(colData(ce)[ce$group == "100 pg RNA, SSIII",] %>% data.frame, aes(TSO, libSizeNormBylib, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~RT_PRIMERS, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_log10( "Normalised counts (arbitrary scale)"
               , breaks = c(0.01, 0.1, 1, 10)) +
  labs( title = "Yield (100 pg RNA, SSIII)") +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_norm_counts_fixed_RNA_facet_RTP-1.png)<!-- -->

Within a library, the reactions are multiplexed by simple pooling without normalisation.  Therefore, reactions that produced more cDNAs will give more sequence reads compared to the others.  Yield increases with TSO amounts, and does not vary much with RT primer amounts, although there may be a trend showing the need to match the RT primer molarity with the TSO molarity.  Effect of RNA amounts and enzyme can not be measured because of the per-library normalisation approach, which confounds with the RNA levels.


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, libSizeNormBylib, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_log10( "Normalised counts (arbitrary scale)"
               , breaks = c(0.01, 0.1, 1, 10)) +
  labs( title = "Sequence yield") +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_norm_counts_by_RTP_facet_RNA-1.png)<!-- -->

### Contour plot


```r
higherIsBetter(z, log(z$libSizeNormBylib), "Relative yield (arbitrary unit on log scale)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/yield_countour-1.png)<!-- -->


## Mapping rate

The percent of proper pairs aligned indicate the amount of data that goes in
the analysis.  The rest is basically discarded.  Here, we do not take tag dust
artefact into account, assuming that they can be removed before sequencing.

The results are somewhat symmetric with the rRNA rate, since rRNA reads are
a large proportion of the discarded data.


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, mapping_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("Mapping rate (% of extracted reads)") +
  ggtitle("Mapping rate")  +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_mapping_rate_by_RTP_facet_RNA-1.png)<!-- -->

### Contour plot


```r
higherIsBetter(z, z$mapping_rate, "Mapping rate (%)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/mapRate_countour-1.png)<!-- -->


## Strand invasion

Strand invasion artefacts are also discarded, but at a later step.  In this
experiment, their amount was reasonably low.

Interestingly, the amount of strand invaders was minimised by high amounts
of TSOs and RT primers.  Does that mean that strand invasion happen first,
and then template-switching happens if primers remain ?


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, strand_invasion_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("Strand invasion rate (% of molecule counts)") +
  labs( title = "Strand invasion"
      , caption = stringr::str_wrap(width = 50, "Strand invasion artefacts are reduced by adding more RT primers.  With SuperScript III, the molarity of TSOs has to be increased at high RNA concentrations.  This does not seem to be teh case with SuperScript IV.")) +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_SI_rate_by_RTP_facet_RNA-1.png)<!-- -->

Because of 1 outlier, the scale is a bit compressed.  Here is the same plot with
the outlier removed.


```r
ggplot(colData(ce[, ce$strand_invasion_rate < 30]) %>% data.frame, aes(TSO, strand_invasion_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("Strand invasion rate (% of molecule counts)") +
  labs( title = "Strand invasion"
      , subtitle = "One outlier removed (1 pg RNA, SSIII, 1 µM RTP, > 30% strand invasion)."
      , caption = stringr::str_wrap(width = 50, "Strand invasion artefacts are reduced by adding more RT primers.  With SuperScript III, the molarity of TSOs has to be increased at high RNA concentrations.  This does not seem to be teh case with SuperScript IV.")) +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_SI_rate_by_RTP_facet_RNA_no_outlier-1.png)<!-- -->

Strand invasion at 100 ng RNA, SSIII, for figure panel

Because of 1 outlier, the scale is a bit compressed.  Here is the same plot with
the outlier removed.


```r
ggplot(subset(colData(ce), group == "100 ng RNA, SSIII") %>% data.frame, aes(TSO, strand_invasion_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("strand invasion (% molecules)") +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_SI_rate_by_RTP_facet_RNA_no_outlier_singlePanel-1.png)<!-- -->

### Contour plot


```r
lowerIsBetter(z, z$strand_invasion_rate, "Strand invasion rate (%)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/SI_countour-1.png)<!-- -->


## Promoter rate

High promoter rate is THE goal of a CAGE experiment.  The molarity of RT
primer does not seem to matter much.  Promoter rate reaches optimum at TSO
molarities higher than 10 µM.

This metric mirrors the strand invasion rate, which suggests that in some
libraries, many strand invaders might be left undetected.



```r
ggplot(colData(ce) %>% data.frame, aes(TSO, promoter_rate, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "fixed", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("Promoter rate (% of molecule counts after SI removal)") +
  labs( title = "Promoter rate"
      , caption = stringr::str_wrap(width = 50, "Promoter rates increases with TSO and RT primer molarity.  It is higher with SuperScript III than with SSIV.")) +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_promoter_rate_by_RTP_facet_RNA-1.png)<!-- -->

Low TSO molarities are much more detrimental for promoter rate at high RNA concentrations.


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, promoter_rate, color=RNA %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~enzyme, scales = "fixed", ncol = 2) +
  scale_color_brewer(name = "RNA (pg)", palette = "YlGnBu") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("Promoter % after SI removal)") +
  labs(title = "Promoter rate") +
  theme_TSO_by_RNA_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_promoter_rate_by_RNA-1.png)<!-- -->

### Contour plots


```r
higherIsBetter(z, z$promoter_rate, "Promoter rate (%)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/promRate_countour-1.png)<!-- -->


```r
zz <- z[z$group.lattice == "10 ng RNA, SSIII",]
higherIsBetter(zz, zz$promoter_rate, "Promoter rate (%)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/promRate_countour_just1panel-1.png)<!-- -->

## Richness (on genes)

Using a free plotting scale as richness reaches extremely low values at 1 pg RNA.

### Richness scale of 10

Many libraries were too shallowly sequenced to allow to calculate richness
on a scale of 100.

Richness is higher when RT primer molarity is higher.


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, r10g, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "free", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TS oligonucleotide molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("Gene richness (on a scale of 10)") +
  labs( title = "Richness on a scale of 10"
      , caption = stringr::str_wrap("Higher RT primer concentration give higher richness.")) +
  theme_TSO_by_RTP_facet_RNA()
```

```
## Warning: Removed 62 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 62 rows containing missing values (geom_point).
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_richness10_by_RTP_facet_RNA-1.png)<!-- -->


```r
higherIsBetter(z, z$r10g, "Richness (on a scale of 10)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/richness10_contour-1.png)<!-- -->

```r
z2 <- z
z2[z2$TSO == 5 & z2$RTP == 4 & z2$group == "1 pg RNA, SSIII", "r10g"] <- NA
higherIsBetter(z2, z2$r10g, "Richness (on a scale of 10)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/richness10_contour-2.png)<!-- -->

### Richness scale of 100

On a scale of 100, we see that TSO concentration needs to be matched with
RNA amounts.


```r
ggplot(colData(ce) %>% data.frame, aes(TSO, r100g, color=RT_PRIMERS %>% factor)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~group, scales = "free", ncol = 2) +
  scale_color_viridis(discrete = TRUE, name = "[RTP] (µM)") +
  scale_x_log10( "TSO molarity (µM)"
               , breaks = c(0.6, 2.5, 10, 40, 160)) +
  scale_y_continuous("Gene richness (on a scale of 100)") +
  labs( title = "Richness on a scale of 100"
      , caption = stringr::str_wrap("Higher RT primer concentration give higher richness.")) +
  theme_TSO_by_RTP_facet_RNA()
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/TSO_vs_richness100_by_RTP_facet_RNA-1.png)<!-- -->


```r
higherIsBetter(z, z$r100g, "Richness (on a scale of 100)")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/richness100_contour-1.png)<!-- -->


Focus on 100 ng RNA, SSIII
==========================


```r
summaryPlot <- function(group) {
  p1 <- lowerIsBetter(z[z$group == group,], z$TD_rate[z$group == group], "TagDust rate (%)", cuts = 5)
    
  p2 <- lowerIsBetter(z[z$group == group,], z$rRNA_rate[z$group == group], "rRNA rate (%)", cuts = 5)
  
  p3 <- higherIsBetter(z[z$group == group,], log(z$libSizeNormBylib[z$group == group]), "Normalised library size (log)", cuts = 5)
   
  p4 <- higherIsBetter(z[z$group == group,], z$mapping_rate[z$group == group], "Mapping rate (%)", cuts = 5)
    
  p5 <- higherIsBetter(z[z$group == group,], z$promoter_rate[z$group == group], "Promoter rate (%)", cuts = 5)
  
  p6 <- higherIsBetter(z[z$group == group,], z$r100g[z$group == group], "Richness (scale 100)", cuts = 5)
    
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
}
summaryPlot("100 ng RNA, SSIII")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/mapRate_contour_lattice-1.png)<!-- -->

```r
summaryPlot("10 pg RNA, SSIV")
```

![](Labcyte-RT_Data_Analysis_8_files/figure-html/mapRate_contour_lattice-2.png)<!-- -->


Session information
===================


```r
sessionInfo()
```

```
## R version 3.5.2 (2018-12-20)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux buster/sid
## 
## Matrix products: default
## BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] viridis_0.5.1               viridisLite_0.3.0          
##  [3] MultiAssayExperiment_1.8.3  SummarizedExperiment_1.12.0
##  [5] DelayedArray_0.8.0          BiocParallel_1.16.6        
##  [7] matrixStats_0.54.0          Biobase_2.42.0             
##  [9] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
## [11] IRanges_2.16.0              S4Vectors_0.20.1           
## [13] BiocGenerics_0.28.0         magrittr_1.5               
## [15] ggplot2_3.0.0               CAGEr_1.24.0               
## 
## loaded via a namespace (and not attached):
##  [1] VGAM_1.1-1               splines_3.5.2           
##  [3] gtools_3.8.1             assertthat_0.2.0        
##  [5] BSgenome_1.50.0          GenomeInfoDbData_1.1.0  
##  [7] Rsamtools_1.34.1         yaml_2.2.0              
##  [9] pillar_1.3.0             lattice_0.20-35         
## [11] glue_1.3.1               digest_0.6.18           
## [13] RColorBrewer_1.1-2       XVector_0.22.0          
## [15] colorspace_1.3-2         htmltools_0.3.6         
## [17] Matrix_1.2-14            plyr_1.8.4              
## [19] XML_3.98-1.16            pkgconfig_2.0.2         
## [21] zlibbioc_1.28.0          purrr_0.2.5             
## [23] scales_1.0.0             stringdist_0.9.5.1      
## [25] tibble_1.4.2             mgcv_1.8-24             
## [27] beanplot_1.2             withr_2.1.2             
## [29] lazyeval_0.2.1           formula.tools_1.7.1     
## [31] crayon_1.3.4             memoise_1.1.0           
## [33] evaluate_0.13            nlme_3.1-137            
## [35] MASS_7.3-50              operator.tools_1.6.3    
## [37] vegan_2.5-4              tools_3.5.2             
## [39] data.table_1.11.4        stringr_1.3.1           
## [41] munsell_0.5.0            cluster_2.0.7-1         
## [43] bindrcpp_0.2.2           Biostrings_2.50.2       
## [45] som_0.3-5.1              compiler_3.5.2          
## [47] rlang_0.3.3              grid_3.5.2              
## [49] RCurl_1.95-4.11          labeling_0.3            
## [51] bitops_1.0-6             rmarkdown_1.12          
## [53] codetools_0.2-15         gtable_0.2.0            
## [55] reshape_0.8.8            R6_2.2.2                
## [57] gridExtra_2.3            GenomicAlignments_1.18.1
## [59] knitr_1.22               dplyr_0.7.6             
## [61] rtracklayer_1.42.2       bindr_0.1.1             
## [63] KernSmooth_2.23-15       permute_0.9-5           
## [65] stringi_1.2.4            Rcpp_0.12.18            
## [67] tidyselect_0.2.4         xfun_0.6
```
