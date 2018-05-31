---
title: "Labcyte-RT Data QC (7th experiment)"
subtitle: "RT optimisation with the Labcyte Echo 525: TSO, RT primer and RNA amounts"
output:
  html_document:
    keep_md: yes
    number_sections: yes
    toc: yes
    pandoc_args:
     - "--lua-filter=links-to-markdown.lua"
---



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

Load summary statistics from MOIRAI
-----------------------------------

MOIRAI shortcuts


```r
MISEQ_RUN      <- "180517_M00528_0364_000000000-BRGK6"
WORKFLOW       <- "OP-WORKFLOW-CAGEscan-short-reads-v2.1.2"
MOIRAI_STAMP   <- "20180521173736"
MOIRAI_PROJ    <- "project/Labcyte_test_decoy" 
MOIRAI_USER    <- "nanoCAGE2" 
ASSEMBLY       <- "mm9"
BASEDIR        <- "/osc-fs_home/scratch/moirai"
MOIRAI_BASE    <- file.path(BASEDIR, MOIRAI_USER)
MOIRAI_RESULTS <- file.path(MOIRAI_BASE, MOIRAI_PROJ, paste(MISEQ_RUN, WORKFLOW, MOIRAI_STAMP, sep = "."))
```


```r
ce <- smallCAGEqc::loadMoiraiStats(
  pipeline  = WORKFLOW,
  multiplex = file.path( MOIRAI_BASE, "input", paste0(MISEQ_RUN, ".multiplex.txt")),
  summary   = file.path( MOIRAI_RESULTS, "text", "summary.txt")) %>% DataFrame

ce$inputFiles <- paste0(MOIRAI_RESULTS, "/CAGEscan_fragments/", ce$samplename, ".bed")

# Discard lines for which input files do not exist.
ce <- ce[sapply(ce$inputFiles, file.exists),]

# Discard lines for which input files are empty.
ce <- ce[file.info(ce$inputFiles)$size != 0,]

ce$inputFilesType <- c("bed")
ce$sampleLabels <- as.character(ce$samplename)

ce <- ce[ce$group != "decoy",]  # Temporary fix

ce
```

```
## DataFrame with 1520 rows and 16 columns
##                      samplename    group  barcode    index     total
##                        <factor> <factor> <factor> <factor> <numeric>
## CTGCGT_TAAGGCGA CTGCGT_TAAGGCGA TAAGGCGA   CTGCGT TAAGGCGA         0
## GCTGCA_TAAGGCGA GCTGCA_TAAGGCGA TAAGGCGA   GCTGCA TAAGGCGA         0
## CGATAC_TAAGGCGA CGATAC_TAAGGCGA TAAGGCGA   CGATAC TAAGGCGA         0
## ACAGAT_TAAGGCGA ACAGAT_TAAGGCGA TAAGGCGA   ACAGAT TAAGGCGA         0
## AGTAGC_TAAGGCGA AGTAGC_TAAGGCGA TAAGGCGA   AGTAGC TAAGGCGA         0
## ...                         ...      ...      ...      ...       ...
## GTAGAT_TCGACGTC GTAGAT_TCGACGTC TCGACGTC   GTAGAT TCGACGTC         0
## AGTCTC_TCGACGTC AGTCTC_TCGACGTC TCGACGTC   AGTCTC TCGACGTC         0
## ATCTAC_TCGACGTC ATCTAC_TCGACGTC TCGACGTC   ATCTAC TCGACGTC         0
## ACATAC_TCGACGTC ACATAC_TCGACGTC TCGACGTC   ACATAC TCGACGTC         0
## CTGCAG_TCGACGTC CTGCAG_TCGACGTC TCGACGTC   CTGCAG TCGACGTC         0
##                 extracted   cleaned   tagdust      rdna    spikes
##                 <numeric> <numeric> <numeric> <numeric> <numeric>
## CTGCGT_TAAGGCGA      2062       783       104      1175         0
## GCTGCA_TAAGGCGA      9540      7596       427      1517         0
## CGATAC_TAAGGCGA      2793      2080       166       544         3
## ACAGAT_TAAGGCGA      1130       541        70       519         0
## AGTAGC_TAAGGCGA     10119      6985       514      2617         3
## ...                   ...       ...       ...       ...       ...
## GTAGAT_TCGACGTC      1120        85      1027         7         1
## AGTCTC_TCGACGTC     35285      2846     32312       124         3
## ATCTAC_TCGACGTC     22731      1584     21098        49         0
## ACATAC_TCGACGTC     24903      2125     22686        90         2
## CTGCAG_TCGACGTC     10539       924      9546        69         0
##                    mapped properpairs    counts
##                 <numeric>   <numeric> <numeric>
## CTGCGT_TAAGGCGA       720         618       603
## GCTGCA_TAAGGCGA      7283        6074      5889
## CGATAC_TAAGGCGA      1982        1618      1571
## ACAGAT_TAAGGCGA       498         384       367
## AGTAGC_TAAGGCGA      6653        5526      5341
## ...                   ...         ...       ...
## GTAGAT_TCGACGTC        13           3         1
## AGTCTC_TCGACGTC       386         183       150
## ATCTAC_TCGACGTC       125          18        16
## ACATAC_TCGACGTC       263         122        88
## CTGCAG_TCGACGTC       178          77        47
##                                                                                                                                                                                                        inputFiles
##                                                                                                                                                                                                       <character>
## CTGCGT_TAAGGCGA /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/CTGCGT_TAAGGCGA.bed
## GCTGCA_TAAGGCGA /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/GCTGCA_TAAGGCGA.bed
## CGATAC_TAAGGCGA /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/CGATAC_TAAGGCGA.bed
## ACAGAT_TAAGGCGA /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/ACAGAT_TAAGGCGA.bed
## AGTAGC_TAAGGCGA /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/AGTAGC_TAAGGCGA.bed
## ...                                                                                                                                                                                                           ...
## GTAGAT_TCGACGTC /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/GTAGAT_TCGACGTC.bed
## AGTCTC_TCGACGTC /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/AGTCTC_TCGACGTC.bed
## ATCTAC_TCGACGTC /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/ATCTAC_TCGACGTC.bed
## ACATAC_TCGACGTC /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/ACATAC_TCGACGTC.bed
## CTGCAG_TCGACGTC /osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte_test_decoy/180517_M00528_0364_000000000-BRGK6.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180521173736/CAGEscan_fragments/CTGCAG_TCGACGTC.bed
##                 inputFilesType    sampleLabels
##                    <character>     <character>
## CTGCGT_TAAGGCGA            bed CTGCGT_TAAGGCGA
## GCTGCA_TAAGGCGA            bed GCTGCA_TAAGGCGA
## CGATAC_TAAGGCGA            bed CGATAC_TAAGGCGA
## ACAGAT_TAAGGCGA            bed ACAGAT_TAAGGCGA
## AGTAGC_TAAGGCGA            bed AGTAGC_TAAGGCGA
## ...                        ...             ...
## GTAGAT_TCGACGTC            bed GTAGAT_TCGACGTC
## AGTCTC_TCGACGTC            bed AGTCTC_TCGACGTC
## ATCTAC_TCGACGTC            bed ATCTAC_TCGACGTC
## ACATAC_TCGACGTC            bed ACATAC_TCGACGTC
## CTGCAG_TCGACGTC            bed CTGCAG_TCGACGTC
```


Load plate design
-----------------

Using transfer designs [6a](Labcyte-RT6a.md), [6b](Labcyte-RT6b.md), [6c](Labcyte-RT6c.md) and [6d](Labcyte-RT6d.md).


```r
plate <- rbind(
  cbind(read.table("plate6a.txt", sep = "\t", header=TRUE, stringsAsFactors=FALSE), repl="Q")
, cbind(read.table("plate6b.txt", sep = "\t", header=TRUE, stringsAsFactors=FALSE), repl="R")
, cbind(read.table("plate6c.txt", sep = "\t", header=TRUE, stringsAsFactors=FALSE), repl="S")
, cbind(read.table("plate6d.txt", sep = "\t", header=TRUE, stringsAsFactors=FALSE), repl="T"))

stopifnot(identical(plate, plate[!duplicated(plate),]))

ce %<>% cbind(plate[match( paste(ce$barcode, ce$index)
                         , paste(plate$BARCODE_SEQ, plate$INDEX)), ])

ce$index %<>% factor(levels = unique(plate$INDEX))  # Keep original order of indexes.
ce$plateID <- ce$repl # Define plate IDs

rownames(ce) %<>% paste(ce$plateID, sep = "_")
ce$sampleLabels <- rownames(ce)

ce$RNA_factor <- ce$RNA %>%
  factor(labels = c("no RNA", "1 pg RNA", "10 pg RNA", "100 pg RNA", "1 ng RNA", "10 ng RNA", "100 ng RNA"))

rm(plate)
```


Create a CAGEexp object and load expression data
------------------------------------------------




```r
if(file.exists(paste0("Labcyte-RT_Data_Analysis_", expNumber, ".Rds"))) {
  ce <- readRDS(paste0("Labcyte-RT_Data_Analysis_", expNumber, ".Rds"))
} else {
  getCTSS(ce, useMulticore = TRUE)
  removeStrandInvaders(ce)
}
```

```
## Loading required namespace: BSgenome.Mmusculus.UCSC.mm9
```


Annotation with GENCODE
-----------------------

Collect annotations and gene symbols via a local GENCODE file
(mm9 GENCODE not available in AnnotationHub)


```r
if(file.exists(paste0("Labcyte-RT_Data_Analysis_", expNumber, ".Rds"))) {
  print("Annotated data loaded from file")
  } else {
  annotateCTSS(ce, rtracklayer::import.gff("/osc-fs_home/scratch/gmtu/annotation/mus_musculus/gencode-M1/gencode.vM1.annotation.gtf.gz"))
  }
```


Save the CAGEexp file
---------------------


```r
saveRDS(ce, paste0("Labcyte-RT_Data_Analysis_", expNumber, ".Rds"))
```


Quality controls
================

Custom _scopes_ displaying _strand invasion_ artefacts.


```r
source("customScopes.R", echo = TRUE)
```

```
## 
## > msScope_qcSI <- function(libs) {
## +     libs$Tag_dust <- libs$extracted - libs$rdna - libs$spikes - 
## +         libs$cleaned
## +     libs$rDNA <- libs$r .... [TRUNCATED] 
## 
## > msScope_counts <- function(libs) {
## +     libs$Promoter <- libs$promoter
## +     libs$Exon <- libs$exon
## +     libs$Intron <- libs$intron
## +     libs$Int .... [TRUNCATED] 
## 
## > msScope_libSizeNormByBarcode <- function(libs) {
## +     libs$Yield <- libs$libSizeNormByBarcode
## +     list(libs = libs, columns = c("Yield"), total = .... [TRUNCATED] 
## 
## > msScope_libSizeNormByIndex <- function(libs) {
## +     libs$Yield <- libs$libSizeNormByIndex
## +     list(libs = libs, columns = c("Yield"), total = lib .... [TRUNCATED] 
## 
## > msScope_libSize <- function(libs) {
## +     libs$Yield <- libs$librarySizes
## +     list(libs = libs, columns = c("Yield"), total = libs$Yield)
## + }
```


Negative controls
-----------------

Some negative controls did not give any sequence and therefore are not in the
CAGEr object at all.


```r
summary(ce$RNA_vol == 0)
```

```
##    Mode   FALSE    TRUE 
## logical    1506      14
```

The signal in negative controls is higher than before in
[experiment 5](Labcyte-RT_Data_Analysis_5.md) and
[experiment 6](Labcyte-RT_Data_Analysis_6.md).  Could that suggest index
swapping?  Not really because in this experiment, a specific set of barcodes
was chosen for negative controls, therefore they are either missed or correct.
Actually, a lot are missed.  However, this is not entirely surprising since
only weak signal is expected.


```r
ce$NC <- ce$RNA_vol == 0

ggpubr::ggarrange( legend = "right", common.legend = TRUE,
  plotAnnot( ce, scope = msScope_qcSI, group = "NC"
           , title = NULL, facet = "index", normalise = FALSE) +
    facet_wrap("facet", ncol = 1) +
    ylab("sequence counts") + xlab("Negative control ?"),
  plotAnnot( ce, scope = msScope_qcSI, group = "NC"
           , title = NULL, facet = "index", normalise = TRUE) +
    facet_wrap("facet", ncol = 1) +
    ylab("Normalised to 100%") + xlab("Negative control ?")
)  %>% ggpubr::annotate_figure(top="QC report, by replicate set")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/neg_ctls_qc-1.png)<!-- -->

The negative controls of `TCCTGAGC`, `GGACTCCT`, `AAGAGGCA`, `GCTCATGA`,
`CGGAGCCT`, `TACGCTGC`m, `ACTGAGCG`, `CCTAAGAC` and `CGATCAGT` have a higher
promoter rate thatn the other samples; this is strange.  It was not the
case in experiment 6.


```r
ggpubr::ggarrange( legend = "right", common.legend = TRUE,
  plotAnnot( ce, scope = msScope_counts, group = "NC"
           , title = NULL, facet = "index", normalise = FALSE) +
    facet_wrap("facet", ncol = 1) +
    ylab("sequence counts") + xlab("Negative control ?"),
  plotAnnot( ce, scope = msScope_counts, group = "NC"
           , title = NULL, facet = "index", normalise = TRUE) +
    facet_wrap("facet", ncol = 1) +
    ylab("Normalised to 100%") + xlab("Negative control ?")
)  %>% ggpubr::annotate_figure(top="QC annotation, by replicate set")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/neg_ctls_annot-1.png)<!-- -->


### Removal of the RNA-negative controls

To ease data handling, the negative controls with no RNA are removed.  By design
all of these negative controls had 10 μM TSO.


```r
colData(ce[,ce$NC]) %>% data.frame %>% summary
```

```
## harmonizing input:
##   removing 1506 sampleMap rows with 'colname' not in colnames of experiments
##   removing 1506 colData rownames not in sampleMap 'primary'
```

```
##            samplename      group      barcode       index       total  
##  ACATCT_TGCAGCTA:1    AAGAGGCA:1   CGATCT :3   TCCTGAGC:1   Min.   :0  
##  ATCAGC_AAGAGGCA:1    ACTGAGCG:1   GTATCT :3   GGACTCCT:1   1st Qu.:0  
##  ATCAGC_CTCTCTAC:1    ATCTCAGG:1   ATCAGC :2   CTCTCTAC:1   Median :0  
##  CACGCA_CCTAAGAC:1    ATGCGCAG:1   CACGCA :2   AAGAGGCA:1   Mean   :0  
##  CACGCA_TACGCTGC:1    CCTAAGAC:1   TCGAGC :2   GCTCATGA:1   3rd Qu.:0  
##  CGATCT_ATCTCAGG:1    CGATCAGT:1   ACATCT :1   ATCTCAGG:1   Max.   :0  
##  (Other)        :8    (Other) :8   (Other):1   (Other) :8              
##    extracted         cleaned         tagdust            rdna       
##  Min.   : 259.0   Min.   : 59.0   Min.   :  12.0   Min.   :  6.00  
##  1st Qu.: 469.8   1st Qu.:138.0   1st Qu.: 138.2   1st Qu.: 37.25  
##  Median : 587.0   Median :178.5   Median : 250.5   Median : 49.50  
##  Mean   : 735.1   Mean   :199.6   Mean   : 425.2   Mean   :109.93  
##  3rd Qu.: 884.0   3rd Qu.:228.2   3rd Qu.: 485.0   3rd Qu.:128.00  
##  Max.   :1890.0   Max.   :380.0   Max.   :1577.0   Max.   :518.00  
##                                                                    
##      spikes           mapped        properpairs         counts      
##  Min.   :0.0000   Min.   : 10.00   Min.   :  4.00   Min.   :  4.00  
##  1st Qu.:0.0000   1st Qu.: 68.75   1st Qu.: 46.75   1st Qu.: 43.25  
##  Median :0.0000   Median :120.50   Median : 92.50   Median : 86.50  
##  Mean   :0.3571   Mean   :135.50   Mean   :106.00   Mean   :101.57  
##  3rd Qu.:0.7500   3rd Qu.:147.25   3rd Qu.:117.75   3rd Qu.:115.50  
##  Max.   :2.0000   Max.   :356.00   Max.   :300.00   Max.   :298.00  
##                                                                     
##   inputFiles        inputFilesType     sampleLabels      
##  Length:14          Length:14          Length:14         
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
##                                                          
##      well               row                 col            sxt           
##  Length:14          Length:14          Min.   : 1.00   Length:14         
##  Class :character   Class :character   1st Qu.:12.00   Class :character  
##  Mode  :character   Mode  :character   Median :15.00   Mode  :character  
##                                        Mean   :14.57                     
##                                        3rd Qu.:19.50                     
##                                        Max.   :24.00                     
##                                                                          
##    BARCODE_ID      TSO_source    TSO_vol        TSO         RT_PRIMERS
##  Min.   : 8.00   Min.   :25   Min.   :25   Min.   :1.25   Min.   :1   
##  1st Qu.:32.00   1st Qu.:25   1st Qu.:25   1st Qu.:1.25   1st Qu.:1   
##  Median :44.00   Median :25   Median :25   Median :1.25   Median :1   
##  Mean   :52.57   Mean   :25   Mean   :25   Mean   :1.25   Mean   :1   
##  3rd Qu.:80.00   3rd Qu.:25   3rd Qu.:25   3rd Qu.:1.25   3rd Qu.:1   
##  Max.   :92.00   Max.   :25   Max.   :25   Max.   :1.25   Max.   :1   
##                                                                       
##  RT_PRIMERS_vol MASTER_MIX_vol    INDEX                RNA       RNA_vol 
##  Min.   :25     Min.   :350    Length:14          Min.   :0   Min.   :0  
##  1st Qu.:25     1st Qu.:350    Class :character   1st Qu.:0   1st Qu.:0  
##  Median :25     Median :350    Mode  :character   Median :0   Median :0  
##  Mean   :25     Mean   :350                       Mean   :0   Mean   :0  
##  3rd Qu.:25     3rd Qu.:350                       3rd Qu.:0   3rd Qu.:0  
##  Max.   :25     Max.   :350                       Max.   :0   Max.   :0  
##                                                                          
##     H2O_vol    BARCODE_SEQ         RNA_level          RTP_level        
##  Min.   :100   Length:14          Length:14          Length:14         
##  1st Qu.:100   Class :character   Class :character   Class :character  
##  Median :100   Mode  :character   Mode  :character   Mode  :character  
##  Mean   :100                                                           
##  3rd Qu.:100                                                           
##  Max.   :100                                                           
##                                                                        
##  plateID repl       RNA_factor  librarySizes    strandInvaders  
##  Q:2     Q:2   no RNA    :14   Min.   :  4.00   Min.   : 0.000  
##  R:4     R:4   1 pg RNA  : 0   1st Qu.: 41.75   1st Qu.: 2.000  
##  S:3     S:3   10 pg RNA : 0   Median : 83.00   Median : 3.000  
##  T:5     T:5   100 pg RNA: 0   Mean   : 97.93   Mean   : 3.643  
##                1 ng RNA  : 0   3rd Qu.:108.50   3rd Qu.: 4.000  
##                10 ng RNA : 0   Max.   :288.00   Max.   :10.000  
##                100 ng RNA: 0                                    
##     promoter           exon            intron          unknown     
##  Min.   :  4.00   Min.   : 0.000   Min.   : 0.000   Min.   : 0.00  
##  1st Qu.: 33.00   1st Qu.: 2.250   1st Qu.: 0.250   1st Qu.: 4.25  
##  Median : 61.00   Median : 4.000   Median : 2.000   Median : 8.50  
##  Mean   : 76.79   Mean   : 5.214   Mean   : 3.143   Mean   :12.79  
##  3rd Qu.: 87.50   3rd Qu.: 8.000   3rd Qu.: 5.250   3rd Qu.:13.75  
##  Max.   :218.00   Max.   :13.000   Max.   :12.000   Max.   :48.00  
##                                                                    
##     NC         
##  Mode:logical  
##  TRUE:14       
##                
##                
##                
##                
## 
```

```r
ce.bak <- ce
ce     <- ce[, ! ce$NC]
```

```
## harmonizing input:
##   removing 14 sampleMap rows with 'colname' not in colnames of experiments
##   removing 14 colData rownames not in sampleMap 'primary'
```


By replicate
------------

Discrepancy between replicates:

 - `CGAGGCTG` has an excess artefacts (tag dust) compared to same libraries from
   other plates.
 - `TAGGCATG` has a different profile compared to other plates, but its amount
   of tag dust clearly indicates that it is correctly identified as a "1 pg"
   library.  Same for `GCTCATGA` in the 10 pg libraries.

The following barcodes were found in libraries where they were not supposed
to be used:

 - `ATCAGC` in `ACTCGCTA`, `GCGTAGTA`, `TAGCGCTC`, `AGGCAGAA`,
`GGAGCTAC`, `GTAGAGGA`, `CGTACTAG` and `TAAGGCGA`.
 - `ACATCT` in `TAGGCATG`.
 - `GTATCT` in `CGAGGCTG`.

In experiment 4, swapped indexes were `CGTACTAG`, `AGGCAGAA` and `TAGGCATG`, 
`CTCTCTAC`.


```r
plotAnnot( ce, scope = msScope_qcSI, group = "RT_PRIMERS", normalise = FALSE
           , title = NULL, facet = "index") + ylab("Number of read pairs") +
  xlab(NULL) + facet_wrap(~facet, ncol = 6, scale = "free")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/qc_plots_by_replicate-1.png)<!-- -->

```r
plotAnnot( ce, scope = msScope_qcSI, group = "RT_PRIMERS", normalise = TRUE
           , title = NULL, facet = "index") + ylab("Number of read pairs") +
  xlab(NULL) + facet_wrap(~facet, ncol = 6, scale = "fixed")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/qc_plots_by_replicate-2.png)<!-- -->

```r
plotAnnot( ce, scope = msScope_qcSI, group = "TSO", normalise = FALSE
           , title = NULL, facet = "index") + ylab("Number of read pairs") +
  xlab(NULL) + facet_wrap(~facet, ncol = 6, scale = "free")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/qc_plots_by_replicate-3.png)<!-- -->


```r
plotAnnot( ce, scope = msScope_qcSI, group = "index", normalise = FALSE
           , title = NULL, facet = "plateID") + ylab("Number of read pairs") +
  xlab(NULL) + facet_wrap(~facet, ncol = 4, scale = "free")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/qc_plots_by_plate-1.png)<!-- -->


By barcode
----------

Barcodes 34 and 7 gave more sequences than average, but their normalised QC
profile or their annotation profile does not stand out


```r
Tier_A <- c(  3, 15, 27, 39, 51, 63, 75
           , 10, 22, 34, 46, 58, 70, 82
           , 11, 23, 35, 47, 59, 71, 83)

Tier_B <- c(  1, 13, 25, 37, 49, 61, 73
           ,  2, 14, 26, 38, 50, 62, 74
           ,  7, 19, 31, 43, 55, 67, 79)

Tier_C <- c(  4, 16, 28, 40, 52, 64, 76
           ,  5, 17, 29, 41, 65, 77, 89
           ,  6, 18, 30, 42, 66, 78, 90)

Tier_N <- c(  8, 20, 32, 44, 68, 80, 92)

ce$BN <- factor(ce$BARCODE_ID, levels = c(Tier_A, Tier_B, Tier_C, Tier_N))

plotAnnot( ce, scope = msScope_qcSI, group = "repl"
           , title = "Sequence counts"
           , facet = "BN", normalise = FALSE) +
  facet_wrap(~facet, ncol=7)
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/plot_by_barcode-1.png)<!-- -->

```r
plotAnnot( ce, scope = msScope_qcSI, group = "repl"
           , title = "Sequence counts (normalised)"
           , facet = "BN", normalise = TRUE) +
  facet_wrap(~facet, ncol=7)
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/plot_by_barcode-2.png)<!-- -->

```r
plotAnnot( ce, scope = msScope_counts, group = "repl"
           , title = "Annotation counts"
           , facet = "BN", normalise = FALSE) +
  facet_wrap(~facet, ncol=7)
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/plot_by_barcode-3.png)<!-- -->


By RT primer quantity
---------------------

Strangely, libraries made with no RT primers have a QC profile that is not
dramatically different from other libraries.  This might again be explained
by contaminations, although the amount of sequences in the "no RT primer"
samples is a bit high for such an explanation.



```r
ggpubr::ggarrange(legend = "right", common.legend = TRUE,
  plotAnnot( ce, scope = msScope_qcSI, group = "RT_PRIMERS", normalise = FALSE
           , title = NULL) + ylab("Number of read pairs") + xlab(NULL),
  plotAnnot( ce, scope = msScope_qcSI, group = "RT_PRIMERS", normalise = TRUE
         , title = NULL) + ylab("Normalised to 100 %") + xlab(NULL)) %>%
  ggpubr::annotate_figure(top="QC control, by amount of RT primers (in μM)")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/qc_rt_primer_quantity-1.png)<!-- -->


```r
ggpubr::ggarrange(legend = "right", common.legend = TRUE,
  plotAnnot( ce, scope = msScope_counts, group = "RT_PRIMERS", normalise = FALSE
           , title = NULL) + ylab("Number of read pairs") + xlab(NULL),
  plotAnnot( ce, scope = msScope_counts, group = "RT_PRIMERS", normalise = TRUE
         , title = NULL) + ylab("Normalised to 100 %") + xlab(NULL)) %>%
  ggpubr::annotate_figure(top="QC annotation, by amount of RT primers (in μM)")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/promrate_rt_primer_quantity-1.png)<!-- -->

Consistency of RTP-negative controls (still looking for index swap)


```r
plotAnnot( ce[,ce$TSO == 10], scope = msScope_qcSI, group = "RT_PRIMERS", normalise = FALSE
           , title = NULL, facet = "index") + ylab("Number of read pairs") +
  xlab(NULL) + facet_wrap(~facet, ncol = 6, scale = "fixed")
```

```
## harmonizing input:
##   removing 1338 sampleMap rows with 'colname' not in colnames of experiments
##   removing 1338 colData rownames not in sampleMap 'primary'
```

```
## Warning: Removed 1344 rows containing missing values (geom_segment).
```

```
## Warning: Removed 1344 rows containing missing values (geom_point).
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

### Removal of the primer-negative controls

To ease data handling (for instance when working with primer ratios), the
negative controls with no primers are removed.


```r
ce <- ce[,ce$RT_PRIMERS != 0]
```

```
## harmonizing input:
##   removing 215 sampleMap rows with 'colname' not in colnames of experiments
##   removing 215 colData rownames not in sampleMap 'primary'
```


## By TSO concentration

 - Highest TSO concentrations give largest molecule counts, except at the
   highest RNA concentration (200 ng / µL)
   
 - Large number of tag artefacts are produced at lowest RNA concentrations
   (20 or 1 pg / µL).  Note that the libraries have not been size-selected
   on purpose, to estimate the amount of artefacts.
   
 - rRNA rate peaks at the highest TSO concentrations, and has another peak
   at intermediate concentrations (~2.5 µM)

 - Strand invation seems to depend more on RNA quantity than on TSO concentration.
 
 -> let's do the next experiment at 10 µM final (200 µM stock), like in the
    nanoCAGE 2017 protocol.  Slightly higher concentrations may have been
    better but we need to consider that our stock plate is soon running out.
    In addition, the actual TSO concentrations in this experiment are lower
    than the expected value, according to the survey of the
    [concentrations of the stock plate PO_8268526](TSO_concentration_check2.md).


```r
ggpubr::ggarrange(legend = "right",
  plotAnnot( ce, scope = msScope_qcSI, group = "TSO", facet="RNA"
           , normalise = FALSE, title = "Yield") +
    ylab("Number of read pairs") + xlab("TSO") +
    facet_wrap(~facet, ncol = 1),
  plotAnnot( ce, scope = msScope_qcSI, group = "TSO", facet="RNA"
           , normalise = TRUE, title = "QC report") +
    ylab("Normalised to 100 %") + xlab(NULL) +
    facet_wrap(~facet, ncol = 1)) %>%
  ggpubr::annotate_figure(top="QC control, by TSO concentration and amount of RNA (ng)")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/qc_TSO_quantity-1.png)<!-- -->

## By RTP concentration

 - At high RNA concentrations, rRNA rate increases with RT primer concentration.
 - At low RNA concentrations, rRNA rate decreases with RT primer concentration...
 
 - at high RNA concentrations, molecule counts only increase marginally with RT
   primer concentration.
   
   -> Let's do the next experiment at 2 µM.


```r
ggpubr::ggarrange(legend = "right",
    plotAnnot( ce, scope = msScope_qcSI, group = "RT_PRIMERS", facet="RNA"
           , normalise = FALSE, title = "Yield") +
    ylab("Number of read pairs") + xlab("TSO") +
    facet_wrap(~facet, ncol = 1),
  plotAnnot( ce, scope = msScope_qcSI, group = "RT_PRIMERS", facet="RNA"
           , normalise = TRUE, title = "QC report") +
    ylab("Normalised to 100 %") + xlab(NULL) +
    facet_wrap(~facet, ncol = 1)) %>%
  ggpubr::annotate_figure(top="QC control, by RTP concentration and amount of RNA (ng)")
```

![](Labcyte-RT_Data_Analysis_7_QC_files/figure-html/qc_norm_TSO_quantity-1.png)<!-- -->


Session information
===================


```r
sessionInfo()
```

```
## R version 3.4.3 (2017-11-30)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## Matrix products: default
## BLAS: /usr/lib/libblas/libblas.so.3.7.0
## LAPACK: /usr/lib/lapack/liblapack.so.3.7.0
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
##  [1] viridis_0.4.0               viridisLite_0.2.0          
##  [3] SummarizedExperiment_1.9.14 DelayedArray_0.4.1         
##  [5] matrixStats_0.52.2          Biobase_2.38.0             
##  [7] GenomicRanges_1.31.19       GenomeInfoDb_1.15.5        
##  [9] IRanges_2.13.26             S4Vectors_0.17.32          
## [11] BiocGenerics_0.25.3         MultiAssayExperiment_1.5.41
## [13] magrittr_1.5                ggplot2_2.2.1              
## [15] CAGEr_1.21.5.1             
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-131                      bitops_1.0-6                     
##  [3] RColorBrewer_1.1-2                rprojroot_1.3-2                  
##  [5] tools_3.4.3                       backports_1.1.2                  
##  [7] R6_2.2.2                          vegan_2.4-5                      
##  [9] platetools_0.0.2                  KernSmooth_2.23-15               
## [11] lazyeval_0.2.1                    mgcv_1.8-22                      
## [13] colorspace_1.3-2                  permute_0.9-4                    
## [15] gridExtra_2.3                     compiler_3.4.3                   
## [17] VennDiagram_1.6.18                labeling_0.3                     
## [19] rtracklayer_1.39.9                scales_0.5.0                     
## [21] stringr_1.3.0                     digest_0.6.15                    
## [23] Rsamtools_1.31.3                  rmarkdown_1.9                    
## [25] stringdist_0.9.4.6                XVector_0.19.8                   
## [27] pkgconfig_2.0.1                   htmltools_0.3.6                  
## [29] BSgenome_1.47.5                   rlang_0.2.0                      
## [31] VGAM_1.0-4                        bindr_0.1                        
## [33] BiocParallel_1.12.0               gtools_3.5.0                     
## [35] dplyr_0.7.4                       RCurl_1.95-4.10                  
## [37] GenomeInfoDbData_0.99.1           futile.logger_1.4.3              
## [39] smallCAGEqc_0.12.2.999999         Matrix_1.2-12                    
## [41] Rcpp_0.12.16                      munsell_0.4.3                    
## [43] stringi_1.1.7                     yaml_2.1.18                      
## [45] MASS_7.3-47                       zlibbioc_1.24.0                  
## [47] plyr_1.8.4                        grid_3.4.3                       
## [49] gdata_2.18.0                      lattice_0.20-35                  
## [51] cowplot_0.9.2                     Biostrings_2.47.9                
## [53] splines_3.4.3                     knitr_1.20                       
## [55] beanplot_1.2                      pillar_1.2.1                     
## [57] ggpubr_0.1.6                      reshape2_1.4.2                   
## [59] codetools_0.2-15                  futile.options_1.0.0             
## [61] XML_3.98-1.9                      glue_1.2.0                       
## [63] evaluate_0.10.1                   lambda.r_1.2                     
## [65] data.table_1.10.4-3               gtable_0.2.0                     
## [67] BSgenome.Mmusculus.UCSC.mm9_1.4.0 purrr_0.2.4                      
## [69] tidyr_0.7.2                       reshape_0.8.7                    
## [71] assertthat_0.2.0                  tibble_1.4.2                     
## [73] som_0.3-5.1                       GenomicAlignments_1.15.12        
## [75] memoise_1.1.0                     bindrcpp_0.2                     
## [77] cluster_2.0.6
```
