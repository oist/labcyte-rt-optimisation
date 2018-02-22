---
title: "Labcyte-RT Data Analysis (merge 3nd experiment)"
output:
  html_document:
    keep_md: yes
    number_sections: yes
    toc: yes
---

Load and merge CAGE data
========================

Data from plates 1~3 and 4~6 were first loaded in separate CAGEexp objects
because our MOIRAI pipeline does not support dual indexing.  Here, we retreive
the objects and merge them.  Merging is relatively slow, so this workflow
just merges the `R` objects and saves the results to a file for further use
by other workflows.


```r
library("CAGEr")
ce1 <- readRDS("Labcyte-RT_Data_Analysis_3a.Rds")
ce2 <- readRDS("Labcyte-RT_Data_Analysis_3b.Rds")
ce  <- mergeCAGEsets(ce1, ce2)
```

```
## Warning in S4Vectors:::set_unlisted_names(unlisted_x, x): failed to set
## names on the unlisted CompressedGRangesList object
```

```r
saveRDS(ce, "Labcyte-RT_Data_Analysis_3.Rds")
```


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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] CAGEr_1.21.5.1
## 
## loaded via a namespace (and not attached):
##  [1] SummarizedExperiment_1.9.14 stringdist_0.9.4.6         
##  [3] gtools_3.5.0                purrr_0.2.4                
##  [5] reshape2_1.4.2              splines_3.4.3              
##  [7] lattice_0.20-35             colorspace_1.3-2           
##  [9] htmltools_0.3.6             stats4_3.4.3               
## [11] mgcv_1.8-22                 beanplot_1.2               
## [13] rtracklayer_1.39.9          yaml_2.1.15                
## [15] XML_3.98-1.9                rlang_0.1.4                
## [17] glue_1.2.0                  BiocParallel_1.12.0        
## [19] BiocGenerics_0.25.3         matrixStats_0.52.2         
## [21] GenomeInfoDbData_0.99.1     plyr_1.8.4                 
## [23] stringr_1.2.0               zlibbioc_1.24.0            
## [25] Biostrings_2.47.9           munsell_0.4.3              
## [27] gtable_0.2.0                VGAM_1.0-4                 
## [29] evaluate_0.10.1             memoise_1.1.0              
## [31] Biobase_2.38.0              knitr_1.17                 
## [33] permute_0.9-4               IRanges_2.13.26            
## [35] MultiAssayExperiment_1.5.41 GenomeInfoDb_1.15.5        
## [37] parallel_3.4.3              Rcpp_0.12.14               
## [39] KernSmooth_2.23-15          backports_1.1.1            
## [41] scales_0.5.0                som_0.3-5.1                
## [43] BSgenome_1.47.5             DelayedArray_0.4.1         
## [45] vegan_2.4-5                 S4Vectors_0.17.32          
## [47] XVector_0.19.8              Rsamtools_1.31.3           
## [49] ggplot2_2.2.1               digest_0.6.12              
## [51] stringi_1.1.6               GenomicRanges_1.31.19      
## [53] grid_3.4.3                  rprojroot_1.2              
## [55] tools_3.4.3                 bitops_1.0-6               
## [57] magrittr_1.5                RCurl_1.95-4.10            
## [59] lazyeval_0.2.1              tibble_1.3.4               
## [61] cluster_2.0.6               tidyr_0.7.2                
## [63] MASS_7.3-47                 Matrix_1.2-12              
## [65] data.table_1.10.4-3         rmarkdown_1.8              
## [67] reshape_0.8.7               nlme_3.1-131               
## [69] GenomicAlignments_1.15.12   compiler_3.4.3
```
