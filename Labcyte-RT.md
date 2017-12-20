---
title: "RT optimisation with the Labcyte Echo 525 (December 2017)"
output: 
  html_document: 
    keep_md: yes
---



Load scripts and libraries.


```r
# We use the R library called "magrittr", to build the data frame describing the
#384-well plate.
library(magrittr)
library(plyr)
# The "setBlock" function to update values in the data frame is shared by multiple
# Rmarkdown files, so it is stored in a separate source file.
source('setBlock.R')
# The "makePlate384" function outputs a data frame representing a 384-well plate
# (one well per row).
source('makePlate384.R')
```

## Plate layout

Create a data frame representing the contents of each well.


```r
plate <- makePlate384()
```

MASTER_MIX
==========

Master mix volume 
------------------

350 nL of mastermix added to each well


```r
plate %<>%
  setBlock("A01~P21", "MASTER_MIX_vol", 350) %>%
  setBlock("K22~L24", "MASTER_MIX_vol", 350) %>%
  setBlock("M23", "MASTER_MIX_vol", 350) %>%
  setBlock("N22~O24", "MASTER_MIX_vol", 350) %>%
  setBlock("P23", "MASTER_MIX_vol", 350)
```

TSO
===

TSO concentration
------------------

 - 80 uM (  +100 nL of 400 uM A1-A7 source wells) in A1-A7, A8-A14, A15-A21, K1-K7, N1-N7 
 - 40 uM ( +50 nL of 400 uM A1-A7 source wells) in B1-B7, B8-B14, B15-B21, L1-L7, O1-O7
 - 20 uM ( +25 nL of 400 uM A1-A7 source wells) in C1-C7, C8-C14, C15-C21, M1-M7, P1-P7
 - 10 uM ( +100 nL of 50 uM source wells) in D1-D7, D8-D14, D15-D21, K8-K14, N8-N14 + J1-J7, J8-J14, J15-J21, K22-K24, L22-L24, M23, N22-N24, O22-O24, P23
 - 5 uM ( +50 nL of 50 uM source wells) in E1-E7, E8-E14, E15-E21, L8-L14, O8-O14
 - 2,5 uM (+25 nL of 50 uM source wells) in F1-F7, F8-F14, F15-F21, M8-M14, P8-P14
 - 1,25 uM (+100 nL of 6,25 uM source wells) in G1-G7, G8-G14, G15-G21, K15-K21, N15-N21
 - 0,625 uM (+50 nL of 6,25 uM source wells) in H1-H7, H8-H14, H15-H21, L15-L21, O15-O21
 - 0,3125 uM (+25 nL of 6,25 uM source wells) in I1-I7, I8-I14, I15-I21, M15-M21, P15-P21
 


```r
plate %<>%
  setBlock("A01~A21", "TSO", 80) %>%
  setBlock("K01~K07", "TSO", 80) %>%
  setBlock("N01~N07", "TSO", 80) %>%
  setBlock("B01~B21", "TSO", 40) %>%
  setBlock("L01~L07", "TSO", 40) %>%
  setBlock("O01~O07", "TSO", 40) %>%
  setBlock("C01~C21", "TSO", 20) %>%
  setBlock("M01~M07", "TSO", 20) %>%
  setBlock("P01~P07", "TSO", 20) %>%
  setBlock("D01~D21", "TSO", 10) %>%
  setBlock("K08~K14", "TSO", 10) %>%
  setBlock("N08~N14", "TSO", 10) %>%
  setBlock("J01~J21", "TSO", 10) %>%
  setBlock("K22~K24", "TSO", 10) %>%
  setBlock("L22~L24", "TSO", 10) %>%
  setBlock("M23", "TSO", 10) %>%
  setBlock("N22~N24", "TSO", 10) %>%
  setBlock("O22~O24", "TSO", 10) %>%
  setBlock("P23", "TSO", 10) %>% 
  setBlock("E01~E21", "TSO", 5) %>%
  setBlock("L08~L14", "TSO", 5) %>%
  setBlock("O08~O14", "TSO", 5) %>%
  setBlock("F01~F21", "TSO", 2.5) %>%
  setBlock("M08~M14", "TSO", 2.5) %>%
  setBlock("P08~P14", "TSO", 2.5) %>%
  setBlock("G01~G21", "TSO", 1.25) %>%
  setBlock("K15~K21", "TSO", 1.25) %>%
  setBlock("N15~N21", "TSO", 1.25) %>%
  setBlock("H01~H21", "TSO", 0.625) %>%
  setBlock("L15~L21", "TSO", 0.625) %>%
  setBlock("O15~O21", "TSO", 0.625) %>%
  setBlock("I01~I21", "TSO", 0.3125) %>%
  setBlock("M15~M21", "TSO", 0.3125) %>%
  setBlock("P15~P21", "TSO", 0.3125)
```

TSO volume
-----------

25 (1-drop), 50 (2-drops) or 100 nL (4-drops) of TSO added, depending on the final concentration


```r
plate %<>%
  setBlock("A01~A21", "TSO_vol", 100) %>%
  setBlock("K01~K07", "TSO_vol", 100) %>%
  setBlock("N01~N07", "TSO_vol", 100) %>%
  setBlock("B01~B21", "TSO_vol", 50) %>%
  setBlock("L01~L07", "TSO_vol", 50) %>%
  setBlock("O01~O07", "TSO_vol", 50) %>%
  setBlock("C01~C21", "TSO_vol", 25) %>%
  setBlock("M01~M07", "TSO_vol", 25) %>%
  setBlock("P01~P07", "TSO_vol", 25) %>%
  setBlock("D01~D21", "TSO_vol", 100) %>%
  setBlock("K08~K14", "TSO_vol", 100) %>%
  setBlock("N08~N14", "TSO_vol", 100) %>%
  setBlock("J01~J21", "TSO_vol", 100) %>%
  setBlock("K22~K24", "TSO_vol", 100) %>%
  setBlock("L22~L24", "TSO_vol", 100) %>%
  setBlock("M23", "TSO_vol", 100) %>%
  setBlock("N22~N24", "TSO_vol", 100) %>%
  setBlock("O22~O24", "TSO_vol", 100) %>%
  setBlock("P23", "TSO_vol", 100) %>% 
  setBlock("E01~E21", "TSO_vol", 50) %>%
  setBlock("L08~L14", "TSO_vol", 50) %>%
  setBlock("O08~O14", "TSO_vol", 50) %>%
  setBlock("F01~F21", "TSO_vol", 25) %>%
  setBlock("M08~M14", "TSO_vol", 25) %>%
  setBlock("P08~P14", "TSO_vol", 25) %>%
  setBlock("G01~G21", "TSO_vol", 100) %>%
  setBlock("K15~K21", "TSO_vol", 100) %>%
  setBlock("N15~N21", "TSO_vol", 100) %>%
  setBlock("H01~H21", "TSO_vol", 50) %>%
  setBlock("L15~L21", "TSO_vol", 50) %>%
  setBlock("O15~O21", "TSO_vol", 50) %>%
  setBlock("I01~I21", "TSO_vol", 25) %>%
  setBlock("M15~M21", "TSO_vol", 25) %>%
  setBlock("P15~P21", "TSO_vol", 25) 
```

Barcode ID
-----------

70 barcodes used for each RNA concentration tested


```r
plate %<>%
  setBlock("A01", "BARCODE_ID", 01) %>%
  setBlock("A02", "BARCODE_ID", 02) %>%
  setBlock("A03", "BARCODE_ID", 03) %>%
  setBlock("A04", "BARCODE_ID", 04) %>%
  setBlock("A05", "BARCODE_ID", 05) %>%
  setBlock("A06", "BARCODE_ID", 06) %>%
  setBlock("A07", "BARCODE_ID", 07) %>%  
  setBlock("A08", "BARCODE_ID", 01) %>%
  setBlock("A09", "BARCODE_ID", 02) %>%
  setBlock("A10", "BARCODE_ID", 03) %>%
  setBlock("A11", "BARCODE_ID", 04) %>%
  setBlock("A12", "BARCODE_ID", 05) %>%
  setBlock("A13", "BARCODE_ID", 06) %>%
  setBlock("A14", "BARCODE_ID", 07) %>%  
  setBlock("A15", "BARCODE_ID", 01) %>%
  setBlock("A16", "BARCODE_ID", 02) %>%
  setBlock("A17", "BARCODE_ID", 03) %>%
  setBlock("A18", "BARCODE_ID", 04) %>%
  setBlock("A19", "BARCODE_ID", 05) %>%
  setBlock("A20", "BARCODE_ID", 06) %>%
  setBlock("A21", "BARCODE_ID", 07) %>%  
  setBlock("K01", "BARCODE_ID", 01) %>%
  setBlock("K02", "BARCODE_ID", 02) %>%
  setBlock("K03", "BARCODE_ID", 03) %>%
  setBlock("K04", "BARCODE_ID", 04) %>%
  setBlock("K05", "BARCODE_ID", 05) %>%
  setBlock("K06", "BARCODE_ID", 06) %>%
  setBlock("K07", "BARCODE_ID", 07) %>%    
  setBlock("N01", "BARCODE_ID", 01) %>%
  setBlock("N02", "BARCODE_ID", 02) %>%
  setBlock("N03", "BARCODE_ID", 03) %>%
  setBlock("N04", "BARCODE_ID", 04) %>%
  setBlock("N05", "BARCODE_ID", 05) %>%
  setBlock("N06", "BARCODE_ID", 06) %>%
  setBlock("N07", "BARCODE_ID", 07) %>%  
  setBlock("B01", "BARCODE_ID", 08) %>%
  setBlock("B02", "BARCODE_ID", 09) %>%
  setBlock("B03", "BARCODE_ID", 10) %>%
  setBlock("B04", "BARCODE_ID", 11) %>%
  setBlock("B05", "BARCODE_ID", 12) %>%
  setBlock("B06", "BARCODE_ID", 13) %>%
  setBlock("B07", "BARCODE_ID", 14) %>%  
  setBlock("B08", "BARCODE_ID", 08) %>%
  setBlock("B09", "BARCODE_ID", 09) %>%
  setBlock("B10", "BARCODE_ID", 10) %>%
  setBlock("B11", "BARCODE_ID", 11) %>%
  setBlock("B12", "BARCODE_ID", 12) %>%
  setBlock("B13", "BARCODE_ID", 13) %>%
  setBlock("B14", "BARCODE_ID", 14) %>%  
  setBlock("B15", "BARCODE_ID", 08) %>%
  setBlock("B16", "BARCODE_ID", 09) %>%
  setBlock("B17", "BARCODE_ID", 10) %>%
  setBlock("B18", "BARCODE_ID", 11) %>%
  setBlock("B19", "BARCODE_ID", 12) %>%
  setBlock("B20", "BARCODE_ID", 13) %>%
  setBlock("B21", "BARCODE_ID", 14) %>% 
  setBlock("L01", "BARCODE_ID", 08) %>%
  setBlock("L02", "BARCODE_ID", 09) %>%
  setBlock("L03", "BARCODE_ID", 10) %>%
  setBlock("L04", "BARCODE_ID", 11) %>%
  setBlock("L05", "BARCODE_ID", 12) %>%
  setBlock("L06", "BARCODE_ID", 13) %>%
  setBlock("L07", "BARCODE_ID", 14) %>% 
  setBlock("O01", "BARCODE_ID", 08) %>%
  setBlock("O02", "BARCODE_ID", 09) %>%
  setBlock("O03", "BARCODE_ID", 10) %>%
  setBlock("O04", "BARCODE_ID", 11) %>%
  setBlock("O05", "BARCODE_ID", 12) %>%
  setBlock("O06", "BARCODE_ID", 13) %>%
  setBlock("O07", "BARCODE_ID", 14) %>% 
  setBlock("C01", "BARCODE_ID", 15) %>%
  setBlock("C02", "BARCODE_ID", 16) %>%
  setBlock("C03", "BARCODE_ID", 17) %>%
  setBlock("C04", "BARCODE_ID", 18) %>%
  setBlock("C05", "BARCODE_ID", 19) %>%
  setBlock("C06", "BARCODE_ID", 20) %>%
  setBlock("C07", "BARCODE_ID", 21) %>%   
  setBlock("C08", "BARCODE_ID", 15) %>%
  setBlock("C09", "BARCODE_ID", 16) %>%
  setBlock("C10", "BARCODE_ID", 17) %>%
  setBlock("C11", "BARCODE_ID", 18) %>%
  setBlock("C12", "BARCODE_ID", 19) %>%
  setBlock("C13", "BARCODE_ID", 20) %>%
  setBlock("C14", "BARCODE_ID", 21) %>%  
  setBlock("C15", "BARCODE_ID", 15) %>%
  setBlock("C16", "BARCODE_ID", 16) %>%
  setBlock("C17", "BARCODE_ID", 17) %>%
  setBlock("C18", "BARCODE_ID", 18) %>%
  setBlock("C19", "BARCODE_ID", 19) %>%
  setBlock("C20", "BARCODE_ID", 20) %>%
  setBlock("C21", "BARCODE_ID", 21) %>%
  setBlock("M01", "BARCODE_ID", 15) %>%
  setBlock("M02", "BARCODE_ID", 16) %>%
  setBlock("M03", "BARCODE_ID", 17) %>%
  setBlock("M04", "BARCODE_ID", 18) %>%
  setBlock("M05", "BARCODE_ID", 19) %>%
  setBlock("M06", "BARCODE_ID", 20) %>%
  setBlock("M07", "BARCODE_ID", 21) %>%
  setBlock("P01", "BARCODE_ID", 15) %>%
  setBlock("P02", "BARCODE_ID", 16) %>%
  setBlock("P03", "BARCODE_ID", 17) %>%
  setBlock("P04", "BARCODE_ID", 18) %>%
  setBlock("P05", "BARCODE_ID", 19) %>%
  setBlock("P06", "BARCODE_ID", 20) %>%
  setBlock("P07", "BARCODE_ID", 21) %>%
  setBlock("D01", "BARCODE_ID", 22) %>%
  setBlock("D02", "BARCODE_ID", 23) %>%
  setBlock("D03", "BARCODE_ID", 24) %>%
  setBlock("D04", "BARCODE_ID", 25) %>%
  setBlock("D05", "BARCODE_ID", 26) %>%
  setBlock("D06", "BARCODE_ID", 27) %>%
  setBlock("D07", "BARCODE_ID", 28) %>%   
  setBlock("D08", "BARCODE_ID", 22) %>%
  setBlock("D09", "BARCODE_ID", 23) %>%
  setBlock("D10", "BARCODE_ID", 24) %>%
  setBlock("D11", "BARCODE_ID", 25) %>%
  setBlock("D12", "BARCODE_ID", 26) %>%
  setBlock("D13", "BARCODE_ID", 27) %>%
  setBlock("D14", "BARCODE_ID", 28) %>% 
  setBlock("D15", "BARCODE_ID", 22) %>%
  setBlock("D16", "BARCODE_ID", 23) %>%
  setBlock("D17", "BARCODE_ID", 24) %>%
  setBlock("D18", "BARCODE_ID", 25) %>%
  setBlock("D19", "BARCODE_ID", 26) %>%
  setBlock("D20", "BARCODE_ID", 27) %>%
  setBlock("D21", "BARCODE_ID", 28) %>%
  setBlock("K08", "BARCODE_ID", 22) %>%
  setBlock("K09", "BARCODE_ID", 23) %>%
  setBlock("K10", "BARCODE_ID", 24) %>%
  setBlock("K11", "BARCODE_ID", 25) %>%
  setBlock("K12", "BARCODE_ID", 26) %>%
  setBlock("K13", "BARCODE_ID", 27) %>%
  setBlock("K14", "BARCODE_ID", 28) %>%  
  setBlock("N08", "BARCODE_ID", 22) %>%
  setBlock("N09", "BARCODE_ID", 23) %>%
  setBlock("N10", "BARCODE_ID", 24) %>%
  setBlock("N11", "BARCODE_ID", 25) %>%
  setBlock("N12", "BARCODE_ID", 26) %>%
  setBlock("N13", "BARCODE_ID", 27) %>%
  setBlock("N14", "BARCODE_ID", 28) %>% 
  setBlock("E01", "BARCODE_ID", 29) %>%
  setBlock("E02", "BARCODE_ID", 30) %>%
  setBlock("E03", "BARCODE_ID", 31) %>%
  setBlock("E04", "BARCODE_ID", 32) %>%
  setBlock("E05", "BARCODE_ID", 33) %>%
  setBlock("E06", "BARCODE_ID", 34) %>%
  setBlock("E07", "BARCODE_ID", 35) %>%
  setBlock("E08", "BARCODE_ID", 29) %>%
  setBlock("E09", "BARCODE_ID", 30) %>%
  setBlock("E10", "BARCODE_ID", 31) %>%
  setBlock("E11", "BARCODE_ID", 32) %>%
  setBlock("E12", "BARCODE_ID", 33) %>%
  setBlock("E13", "BARCODE_ID", 34) %>%
  setBlock("E14", "BARCODE_ID", 35) %>%  
  setBlock("E15", "BARCODE_ID", 29) %>%
  setBlock("E16", "BARCODE_ID", 30) %>%
  setBlock("E17", "BARCODE_ID", 31) %>%
  setBlock("E18", "BARCODE_ID", 32) %>%
  setBlock("E19", "BARCODE_ID", 33) %>%
  setBlock("E20", "BARCODE_ID", 34) %>%
  setBlock("E21", "BARCODE_ID", 35) %>%  
  setBlock("L08", "BARCODE_ID", 29) %>%
  setBlock("L09", "BARCODE_ID", 30) %>%
  setBlock("L10", "BARCODE_ID", 31) %>%
  setBlock("L11", "BARCODE_ID", 32) %>%
  setBlock("L12", "BARCODE_ID", 33) %>%
  setBlock("L13", "BARCODE_ID", 34) %>%
  setBlock("L14", "BARCODE_ID", 35) %>%  
  setBlock("O08", "BARCODE_ID", 29) %>%
  setBlock("O09", "BARCODE_ID", 30) %>%
  setBlock("O10", "BARCODE_ID", 31) %>%
  setBlock("O11", "BARCODE_ID", 32) %>%
  setBlock("O12", "BARCODE_ID", 33) %>%
  setBlock("O13", "BARCODE_ID", 34) %>%
  setBlock("O14", "BARCODE_ID", 35) %>%    
  setBlock("F01", "BARCODE_ID", 36) %>%
  setBlock("F02", "BARCODE_ID", 37) %>%
  setBlock("F03", "BARCODE_ID", 38) %>%
  setBlock("F04", "BARCODE_ID", 39) %>%
  setBlock("F05", "BARCODE_ID", 40) %>%
  setBlock("F06", "BARCODE_ID", 41) %>%
  setBlock("F07", "BARCODE_ID", 42) %>% 
  setBlock("F08", "BARCODE_ID", 36) %>%
  setBlock("F09", "BARCODE_ID", 37) %>%
  setBlock("F10", "BARCODE_ID", 38) %>%
  setBlock("F11", "BARCODE_ID", 39) %>%
  setBlock("F12", "BARCODE_ID", 40) %>%
  setBlock("F13", "BARCODE_ID", 41) %>%
  setBlock("F14", "BARCODE_ID", 42) %>%  
  setBlock("F15", "BARCODE_ID", 36) %>%
  setBlock("F16", "BARCODE_ID", 37) %>%
  setBlock("F17", "BARCODE_ID", 38) %>%
  setBlock("F18", "BARCODE_ID", 39) %>%
  setBlock("F19", "BARCODE_ID", 40) %>%
  setBlock("F20", "BARCODE_ID", 41) %>%
  setBlock("F21", "BARCODE_ID", 42) %>%
  setBlock("M08", "BARCODE_ID", 36) %>%
  setBlock("M09", "BARCODE_ID", 37) %>%
  setBlock("M10", "BARCODE_ID", 38) %>%
  setBlock("M11", "BARCODE_ID", 39) %>%
  setBlock("M12", "BARCODE_ID", 40) %>%
  setBlock("M13", "BARCODE_ID", 41) %>%
  setBlock("M14", "BARCODE_ID", 42) %>%   
  setBlock("P08", "BARCODE_ID", 36) %>%
  setBlock("P09", "BARCODE_ID", 37) %>%
  setBlock("P10", "BARCODE_ID", 38) %>%
  setBlock("P11", "BARCODE_ID", 39) %>%
  setBlock("P12", "BARCODE_ID", 40) %>%
  setBlock("P13", "BARCODE_ID", 41) %>%
  setBlock("P14", "BARCODE_ID", 42) %>%  
  setBlock("G01", "BARCODE_ID", 43) %>%
  setBlock("G02", "BARCODE_ID", 44) %>%
  setBlock("G03", "BARCODE_ID", 45) %>%
  setBlock("G04", "BARCODE_ID", 46) %>%
  setBlock("G05", "BARCODE_ID", 47) %>%
  setBlock("G06", "BARCODE_ID", 48) %>%
  setBlock("G07", "BARCODE_ID", 49) %>% 
  setBlock("G08", "BARCODE_ID", 43) %>%
  setBlock("G09", "BARCODE_ID", 44) %>%
  setBlock("G10", "BARCODE_ID", 45) %>%
  setBlock("G11", "BARCODE_ID", 46) %>%
  setBlock("G12", "BARCODE_ID", 47) %>%
  setBlock("G13", "BARCODE_ID", 48) %>%
  setBlock("G14", "BARCODE_ID", 49) %>%   
  setBlock("G15", "BARCODE_ID", 43) %>%
  setBlock("G16", "BARCODE_ID", 44) %>%
  setBlock("G17", "BARCODE_ID", 45) %>%
  setBlock("G18", "BARCODE_ID", 46) %>%
  setBlock("G19", "BARCODE_ID", 47) %>%
  setBlock("G20", "BARCODE_ID", 48) %>%
  setBlock("G21", "BARCODE_ID", 49) %>% 
  setBlock("K15", "BARCODE_ID", 43) %>%
  setBlock("K16", "BARCODE_ID", 44) %>%
  setBlock("K17", "BARCODE_ID", 45) %>%
  setBlock("K18", "BARCODE_ID", 46) %>%
  setBlock("K19", "BARCODE_ID", 47) %>%
  setBlock("K20", "BARCODE_ID", 48) %>%
  setBlock("K21", "BARCODE_ID", 49) %>%   
  setBlock("N15", "BARCODE_ID", 43) %>%
  setBlock("N16", "BARCODE_ID", 44) %>%
  setBlock("N17", "BARCODE_ID", 45) %>%
  setBlock("N18", "BARCODE_ID", 46) %>%
  setBlock("N19", "BARCODE_ID", 47) %>%
  setBlock("N20", "BARCODE_ID", 48) %>%
  setBlock("N21", "BARCODE_ID", 49) %>%  
  setBlock("H01", "BARCODE_ID", 50) %>%
  setBlock("H02", "BARCODE_ID", 51) %>%
  setBlock("H03", "BARCODE_ID", 52) %>%
  setBlock("H04", "BARCODE_ID", 53) %>%
  setBlock("H05", "BARCODE_ID", 54) %>%
  setBlock("H06", "BARCODE_ID", 55) %>%
  setBlock("H07", "BARCODE_ID", 56) %>%
  setBlock("H08", "BARCODE_ID", 50) %>%
  setBlock("H09", "BARCODE_ID", 51) %>%
  setBlock("H10", "BARCODE_ID", 52) %>%
  setBlock("H11", "BARCODE_ID", 53) %>%
  setBlock("H12", "BARCODE_ID", 54) %>%
  setBlock("H13", "BARCODE_ID", 55) %>%
  setBlock("H14", "BARCODE_ID", 56) %>%
  setBlock("H15", "BARCODE_ID", 50) %>%
  setBlock("H16", "BARCODE_ID", 51) %>%
  setBlock("H17", "BARCODE_ID", 52) %>%
  setBlock("H18", "BARCODE_ID", 53) %>%
  setBlock("H19", "BARCODE_ID", 54) %>%
  setBlock("H20", "BARCODE_ID", 55) %>%
  setBlock("H21", "BARCODE_ID", 56) %>%
  setBlock("L15", "BARCODE_ID", 50) %>%
  setBlock("L16", "BARCODE_ID", 51) %>%
  setBlock("L17", "BARCODE_ID", 52) %>%
  setBlock("L18", "BARCODE_ID", 53) %>%
  setBlock("L19", "BARCODE_ID", 54) %>%
  setBlock("L20", "BARCODE_ID", 55) %>%
  setBlock("L21", "BARCODE_ID", 56) %>%  
  setBlock("O15", "BARCODE_ID", 50) %>%
  setBlock("O16", "BARCODE_ID", 51) %>%
  setBlock("O17", "BARCODE_ID", 52) %>%
  setBlock("O18", "BARCODE_ID", 53) %>%
  setBlock("O19", "BARCODE_ID", 54) %>%
  setBlock("O20", "BARCODE_ID", 55) %>%
  setBlock("O21", "BARCODE_ID", 56) %>% 
  setBlock("I01", "BARCODE_ID", 57) %>%
  setBlock("I02", "BARCODE_ID", 58) %>%
  setBlock("I03", "BARCODE_ID", 59) %>%
  setBlock("I04", "BARCODE_ID", 60) %>%
  setBlock("I05", "BARCODE_ID", 61) %>%
  setBlock("I06", "BARCODE_ID", 62) %>%
  setBlock("I07", "BARCODE_ID", 63) %>%  
  setBlock("I08", "BARCODE_ID", 57) %>%
  setBlock("I09", "BARCODE_ID", 58) %>%
  setBlock("I10", "BARCODE_ID", 59) %>%
  setBlock("I11", "BARCODE_ID", 60) %>%
  setBlock("I12", "BARCODE_ID", 61) %>%
  setBlock("I13", "BARCODE_ID", 62) %>%
  setBlock("I14", "BARCODE_ID", 63) %>%    
  setBlock("I15", "BARCODE_ID", 57) %>%
  setBlock("I16", "BARCODE_ID", 58) %>%
  setBlock("I17", "BARCODE_ID", 59) %>%
  setBlock("I18", "BARCODE_ID", 60) %>%
  setBlock("I19", "BARCODE_ID", 61) %>%
  setBlock("I20", "BARCODE_ID", 62) %>%
  setBlock("I21", "BARCODE_ID", 63) %>%   
  setBlock("M15", "BARCODE_ID", 57) %>%
  setBlock("M16", "BARCODE_ID", 58) %>%
  setBlock("M17", "BARCODE_ID", 59) %>%
  setBlock("M18", "BARCODE_ID", 60) %>%
  setBlock("M19", "BARCODE_ID", 61) %>%
  setBlock("M20", "BARCODE_ID", 62) %>%
  setBlock("M21", "BARCODE_ID", 63) %>%  
  setBlock("P15", "BARCODE_ID", 57) %>%
  setBlock("P16", "BARCODE_ID", 58) %>%
  setBlock("P17", "BARCODE_ID", 59) %>%
  setBlock("P18", "BARCODE_ID", 60) %>%
  setBlock("P19", "BARCODE_ID", 61) %>%
  setBlock("P20", "BARCODE_ID", 62) %>%
  setBlock("P21", "BARCODE_ID", 63) %>%   
  setBlock("J01", "BARCODE_ID", 64) %>%
  setBlock("J02", "BARCODE_ID", 65) %>%
  setBlock("J03", "BARCODE_ID", 66) %>%
  setBlock("J04", "BARCODE_ID", 67) %>%
  setBlock("J05", "BARCODE_ID", 68) %>%
  setBlock("J06", "BARCODE_ID", 69) %>%
  setBlock("J07", "BARCODE_ID", 70) %>%   
  setBlock("J08", "BARCODE_ID", 64) %>%
  setBlock("J09", "BARCODE_ID", 65) %>%
  setBlock("J10", "BARCODE_ID", 66) %>%
  setBlock("J11", "BARCODE_ID", 67) %>%
  setBlock("J12", "BARCODE_ID", 68) %>%
  setBlock("J13", "BARCODE_ID", 69) %>%
  setBlock("J14", "BARCODE_ID", 70) %>%   
  setBlock("J15", "BARCODE_ID", 64) %>%
  setBlock("J16", "BARCODE_ID", 65) %>%
  setBlock("J17", "BARCODE_ID", 66) %>%
  setBlock("J18", "BARCODE_ID", 67) %>%
  setBlock("J19", "BARCODE_ID", 68) %>%
  setBlock("J20", "BARCODE_ID", 69) %>%
  setBlock("J21", "BARCODE_ID", 70) %>%  
  setBlock("M23", "BARCODE_ID", 64) %>%
  setBlock("K22", "BARCODE_ID", 65) %>%
  setBlock("K23", "BARCODE_ID", 66) %>%
  setBlock("K24", "BARCODE_ID", 67) %>%
  setBlock("L22", "BARCODE_ID", 68) %>%
  setBlock("L23", "BARCODE_ID", 69) %>%
  setBlock("L24", "BARCODE_ID", 70) %>%  
  setBlock("P23", "BARCODE_ID", 64) %>%
  setBlock("N22", "BARCODE_ID", 65) %>%
  setBlock("N23", "BARCODE_ID", 66) %>%
  setBlock("N24", "BARCODE_ID", 67) %>%
  setBlock("O22", "BARCODE_ID", 68) %>%
  setBlock("O23", "BARCODE_ID", 69) %>%
  setBlock("O24", "BARCODE_ID", 70)   
```

RT_PRIMERS
===========

RT primers concentration
-------------------------

 - 0 uM (controls, + 0 nL) in Col 1, 8 and 15 + M23 + P23
 - 0.125 uM (+ 25 nL from source well K1 2.5 uM) in Col 2, 9 and 16 + K22 + N22
 - 0.25 uM (+ 25 nL from source well K2 5 uM) in Col 3, 10 and 17 + K23 + N23
 - 0.5 uM (+ 25 nL from source well K3 10 uM) in Col 4, 11 and 18 + K24 + N24
 - 1 uM (+ 25 nL from source well K4 20 uM) in Col 5, 12 and 19 + L22 + O22
 - 2 uM (+ 25 nL from source well K5 40 uM) in Col 6, 13 and 20 + L23 + O23
 - 4 uM (+ 25 nL from source well K6 80 uM) in Col 7, 14 and 21 + L24 + O24


```r
plate %<>%
  setBlock("A01~P01",   "RT_PRIMERS", 0) %>%
  setBlock("A08~P08",   "RT_PRIMERS", 0) %>%
  setBlock("A15~P15",   "RT_PRIMERS", 0) %>%
  setBlock("M23",   "RT_PRIMERS", 0) %>%
  setBlock("P23",   "RT_PRIMERS", 0) %>%
  setBlock("A02~P02",   "RT_PRIMERS", 0.125) %>%
  setBlock("A09~P09",   "RT_PRIMERS", 0.125) %>%
  setBlock("A16~P16",   "RT_PRIMERS", 0.125) %>%
  setBlock("K22",   "RT_PRIMERS", 0.125) %>%
  setBlock("N22",   "RT_PRIMERS", 0.125) %>%
  setBlock("A03~P03",   "RT_PRIMERS", 0.25) %>%
  setBlock("A10~P10",   "RT_PRIMERS", 0.25) %>%
  setBlock("A17~P17",   "RT_PRIMERS", 0.25) %>%
  setBlock("K23",   "RT_PRIMERS", 0.25) %>%
  setBlock("N23",   "RT_PRIMERS", 0.25) %>%
  setBlock("A04~P04",   "RT_PRIMERS", 0.5) %>%
  setBlock("A11~P11",   "RT_PRIMERS", 0.5) %>%
  setBlock("A18~P18",   "RT_PRIMERS", 0.5) %>%
  setBlock("K24",   "RT_PRIMERS", 0.5) %>%
  setBlock("N24",   "RT_PRIMERS", 0.5) %>%
  setBlock("A05~P05",   "RT_PRIMERS", 1) %>%
  setBlock("A12~P12",   "RT_PRIMERS", 1) %>%
  setBlock("A19~P19",   "RT_PRIMERS", 1) %>%
  setBlock("L22",   "RT_PRIMERS", 1) %>%
  setBlock("O22",   "RT_PRIMERS", 1) %>%
  setBlock("A06~P06",   "RT_PRIMERS", 2) %>%
  setBlock("A13~P13",   "RT_PRIMERS", 2) %>%
  setBlock("A20~P20",   "RT_PRIMERS", 2) %>%
  setBlock("L23",   "RT_PRIMERS", 2) %>%
  setBlock("O23",   "RT_PRIMERS", 2) %>%
  setBlock("A07~P07",   "RT_PRIMERS", 4) %>%
  setBlock("A14~P14",   "RT_PRIMERS", 4) %>%
  setBlock("A21~P21",   "RT_PRIMERS", 4) %>%
  setBlock("L24",   "RT_PRIMERS", 4) %>%
  setBlock("O24",   "RT_PRIMERS", 4) 
```

RT primer volume 
------------------

25 nl of RT_PRIMERS added in each well, except the negative controls 


```r
plate %<>%
  setBlock("A01~P01",   "RT_PRIMERS_vol", 0) %>%
  setBlock("A08~P08",   "RT_PRIMERS_vol", 0) %>%
  setBlock("A15~P15",   "RT_PRIMERS_vol", 0) %>%
  setBlock("M23",   "RT_PRIMERS_vol", 0) %>%
  setBlock("P23",   "RT_PRIMERS_vol", 0) %>%
  setBlock("A02~P02",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A09~P09",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A16~P16",   "RT_PRIMERS_vol", 25) %>%
  setBlock("K22",   "RT_PRIMERS_vol", 25) %>%
  setBlock("N22",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A03~P03",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A10~P10",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A17~P17",   "RT_PRIMERS_vol", 25) %>%
  setBlock("K23",   "RT_PRIMERS_vol", 25) %>%
  setBlock("N23",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A04~P04",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A11~P11",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A18~P18",   "RT_PRIMERS_vol", 25) %>%
  setBlock("K24",   "RT_PRIMERS_vol", 25) %>%
  setBlock("N24",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A05~P05",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A12~P12",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A19~P19",   "RT_PRIMERS_vol", 25) %>%
  setBlock("L22",   "RT_PRIMERS_vol", 25) %>%
  setBlock("O22",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A06~P06",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A13~P13",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A20~P20",   "RT_PRIMERS_vol", 25) %>%
  setBlock("L23",   "RT_PRIMERS_vol", 25) %>%
  setBlock("O23",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A07~P07",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A14~P14",   "RT_PRIMERS_vol", 25) %>%
  setBlock("A21~P21",   "RT_PRIMERS_vol", 25) %>%
  setBlock("L24",   "RT_PRIMERS_vol", 25) %>%
  setBlock("O24",   "RT_PRIMERS_vol", 25) 
```

RNA
====

RNA concentration
------------------

 - 0 (controls, + 0 nL) in J01-J21 + K22-K24 + L22-L24 + M23 + N22-N24 + O22-O24 + P23
 - 100 ng (+ 25 nL from source well L1 4 ug/uL) in A01-I07
 - 10 ng (+ 25 nL from source well L2 400 ng/uL) in A08-I14
 - 1 ng (+ 25 nL from source well L3 40 ng/uL) in A15-I21
 - 100 pg (+ 25 nL from source well L4 4 ng/uL) in K01-M21
 - 10 pg (+ 25 nL from source well L5 0.4 ng/uL) in N01-P21


```r
plate %<>%
  setBlock("J01~J21", "RNA", 0) %>%
  setBlock("K22~L24", "RNA", 0) %>%
  setBlock("M23", "RNA", 0) %>%
  setBlock("N22~O24", "RNA", 0) %>%
  setBlock("P23", "RNA", 0) %>%
  setBlock("A01~I07", "RNA", 100) %>%
  setBlock("A08~I14", "RNA", 10) %>%
  setBlock("A15~I21", "RNA", 1) %>%
  setBlock("K01~M21", "RNA", 0.1) %>%
  setBlock("N01~P21", "RNA", 0.01)
```

RNA volume
-----------

25 nL of RNA added to each well, except the negative controls


```r
plate %<>%
  setBlock("J01~J21", "RNA_vol", 0) %>%
  setBlock("K22~L24", "RNA_vol", 0) %>%
  setBlock("M23", "RNA_vol", 0) %>%
  setBlock("N22~O24", "RNA_vol", 0) %>%
  setBlock("P23", "RNA_vol", 0) %>%
  setBlock("A01~I07", "RNA_vol", 25) %>%
  setBlock("A08~I14", "RNA_vol", 25) %>%
  setBlock("A15~I21", "RNA_vol", 25) %>%
  setBlock("K01~M21", "RNA_vol", 25) %>%
  setBlock("N01~P21", "RNA_vol", 25)
```

H2O
====

H2O volume
-----------

0, 25, 50, 75 or 100 nL H2O added to complete RT reaction volume up to 500 nL


```r
plate %<>%
  setBlock("A02~A07", "H2O_vol", 0) %>%
  setBlock("A09~A14", "H2O_vol", 0) %>%
  setBlock("A16~A21", "H2O_vol", 0) %>%
  setBlock("D02~D07", "H2O_vol", 0) %>%
  setBlock("D09~D14", "H2O_vol", 0) %>%
  setBlock("D16~D21", "H2O_vol", 0) %>%
  setBlock("G02~G07", "H2O_vol", 0) %>%
  setBlock("G09~G14", "H2O_vol", 0) %>%
  setBlock("G16~G21", "H2O_vol", 0) %>%
  setBlock("K02~K07", "H2O_vol", 0) %>%
  setBlock("K09~K14", "H2O_vol", 0) %>%
  setBlock("K16~K21", "H2O_vol", 0) %>%
  setBlock("N02~N07", "H2O_vol", 0) %>%
  setBlock("N09~N14", "H2O_vol", 0) %>%
  setBlock("N16~N21", "H2O_vol", 0) %>%
  setBlock("J02~J07", "H2O_vol", 25) %>%
  setBlock("J09~J14", "H2O_vol", 25) %>%
  setBlock("J16~J21", "H2O_vol", 25) %>%
  setBlock("J02~J07", "H2O_vol", 25) %>%
  setBlock("K22~K24", "H2O_vol", 25) %>%
  setBlock("L22~L24", "H2O_vol", 25) %>%
  setBlock("N22~N24", "H2O_vol", 25) %>%
  setBlock("O22~O24", "H2O_vol", 25) %>%
  setBlock("K22~K24", "H2O_vol", 25) %>%
  setBlock("A01", "H2O_vol", 25) %>%
  setBlock("A08", "H2O_vol", 25) %>%
  setBlock("A15", "H2O_vol", 25) %>%
  setBlock("D01", "H2O_vol", 25) %>%
  setBlock("D08", "H2O_vol", 25) %>%
  setBlock("D15", "H2O_vol", 25) %>%
  setBlock("G01", "H2O_vol", 25) %>%
  setBlock("G08", "H2O_vol", 25) %>%
  setBlock("G15", "H2O_vol", 25) %>%
  setBlock("K01", "H2O_vol", 25) %>%
  setBlock("K08", "H2O_vol", 25) %>%
  setBlock("K15", "H2O_vol", 25) %>%
  setBlock("N01", "H2O_vol", 25) %>%
  setBlock("N08", "H2O_vol", 25) %>%
  setBlock("N15", "H2O_vol", 25) %>%
  setBlock("B02~B07", "H2O_vol", 50) %>%
  setBlock("B09~B14", "H2O_vol", 50) %>%
  setBlock("B16~B21", "H2O_vol", 50) %>%
  setBlock("E02~E07", "H2O_vol", 50) %>%
  setBlock("E09~E14", "H2O_vol", 50) %>%
  setBlock("E16~E21", "H2O_vol", 50) %>%
  setBlock("H02~H07", "H2O_vol", 50) %>%
  setBlock("H09~H14", "H2O_vol", 50) %>%
  setBlock("H16~H21", "H2O_vol", 50) %>%
  setBlock("L02~L07", "H2O_vol", 50) %>%
  setBlock("L09~L14", "H2O_vol", 50) %>%
  setBlock("L16~L21", "H2O_vol", 50) %>%
  setBlock("O02~O07", "H2O_vol", 50) %>%
  setBlock("O09~O14", "H2O_vol", 50) %>%
  setBlock("O16~O21", "H2O_vol", 50) %>%
  setBlock("M23", "H2O_vol", 50) %>%
  setBlock("P23", "H2O_vol", 50) %>%
  setBlock("J01", "H2O_vol", 50) %>%
  setBlock("J08", "H2O_vol", 50) %>%
  setBlock("J15", "H2O_vol", 50) %>%
  setBlock("C02~C07", "H2O_vol", 75) %>%
  setBlock("C09~C14", "H2O_vol", 75) %>%
  setBlock("C16~C21", "H2O_vol", 75) %>%
  setBlock("F02~F07", "H2O_vol", 75) %>%
  setBlock("F09~F14", "H2O_vol", 75) %>%
  setBlock("F16~F21", "H2O_vol", 75) %>%
  setBlock("I02~I07", "H2O_vol", 75) %>%
  setBlock("I09~I14", "H2O_vol", 75) %>%
  setBlock("I16~I21", "H2O_vol", 75) %>%
  setBlock("M02~M07", "H2O_vol", 75) %>%
  setBlock("M09~M14", "H2O_vol", 75) %>%
  setBlock("M16~M21", "H2O_vol", 75) %>%
  setBlock("P02~P07", "H2O_vol", 75) %>%
  setBlock("P09~P14", "H2O_vol", 75) %>%
  setBlock("P16~P21", "H2O_vol", 75) %>%
  setBlock("B01", "H2O_vol", 75) %>%
  setBlock("B08", "H2O_vol", 75) %>%
  setBlock("B15", "H2O_vol", 75) %>%
  setBlock("E01", "H2O_vol", 75) %>%
  setBlock("E08", "H2O_vol", 75) %>%
  setBlock("E15", "H2O_vol", 75) %>%
  setBlock("H01", "H2O_vol", 75) %>%
  setBlock("H08", "H2O_vol", 75) %>%
  setBlock("H15", "H2O_vol", 75) %>%
  setBlock("L01", "H2O_vol", 75) %>%
  setBlock("L08", "H2O_vol", 75) %>%
  setBlock("L15", "H2O_vol", 75) %>%
  setBlock("O01", "H2O_vol", 75) %>%
  setBlock("O08", "H2O_vol", 75) %>%
  setBlock("O15", "H2O_vol", 75) %>%
  setBlock("C01", "H2O_vol", 100) %>%
  setBlock("C08", "H2O_vol", 100) %>%
  setBlock("C15", "H2O_vol", 100) %>%
  setBlock("F01", "H2O_vol", 100) %>%
  setBlock("F08", "H2O_vol", 100) %>%
  setBlock("F15", "H2O_vol", 100) %>%
  setBlock("I01", "H2O_vol", 100) %>%
  setBlock("I08", "H2O_vol", 100) %>%
  setBlock("I15", "H2O_vol", 100) %>%
  setBlock("M01", "H2O_vol", 100) %>%
  setBlock("M08", "H2O_vol", 100) %>%
  setBlock("M15", "H2O_vol", 100) %>%
  setBlock("P01", "H2O_vol", 100) %>%
  setBlock("P08", "H2O_vol", 100) %>%
  setBlock("P15", "H2O_vol", 100)
```

RATIO TSO/RT_PRIMERS
=====================

Different ratio of TSO/RT_PRIMERS tested


```r
plate$PRIMERS_RATIO <- ""
plate$PRIMERS_RATIO <- c(plate$TSO/plate$RT_PRIMERS)
plate$PRIMERS_RATIO <- sub("Inf", "no_RT_PRIMERS", plate$PRIMERS_RATIO)
```


```r
plate$total_volume <- ""
plate$total_volume <- rowSums(plate[, c(3,5,8,10,11)])
plate 
```

```
##     Row Col MASTER_MIX_vol     TSO TSO_vol BARCODE_ID RT_PRIMERS RT_PRIMERS_vol   RNA RNA_vol H2O_vol PRIMERS_RATIO
## A1    A   1            350 80.0000     100          1      0.000              0 1e+02      25      25 no_RT_PRIMERS
## B1    B   1            350 40.0000      50          8      0.000              0 1e+02      25      75 no_RT_PRIMERS
## C1    C   1            350 20.0000      25         15      0.000              0 1e+02      25     100 no_RT_PRIMERS
## D1    D   1            350 10.0000     100         22      0.000              0 1e+02      25      25 no_RT_PRIMERS
## E1    E   1            350  5.0000      50         29      0.000              0 1e+02      25      75 no_RT_PRIMERS
## F1    F   1            350  2.5000      25         36      0.000              0 1e+02      25     100 no_RT_PRIMERS
## G1    G   1            350  1.2500     100         43      0.000              0 1e+02      25      25 no_RT_PRIMERS
## H1    H   1            350  0.6250      50         50      0.000              0 1e+02      25      75 no_RT_PRIMERS
## I1    I   1            350  0.3125      25         57      0.000              0 1e+02      25     100 no_RT_PRIMERS
## J1    J   1            350 10.0000     100         64      0.000              0 0e+00       0      50 no_RT_PRIMERS
## K1    K   1            350 80.0000     100          1      0.000              0 1e-01      25      25 no_RT_PRIMERS
## L1    L   1            350 40.0000      50          8      0.000              0 1e-01      25      75 no_RT_PRIMERS
## M1    M   1            350 20.0000      25         15      0.000              0 1e-01      25     100 no_RT_PRIMERS
## N1    N   1            350 80.0000     100          1      0.000              0 1e-02      25      25 no_RT_PRIMERS
## O1    O   1            350 40.0000      50          8      0.000              0 1e-02      25      75 no_RT_PRIMERS
## P1    P   1            350 20.0000      25         15      0.000              0 1e-02      25     100 no_RT_PRIMERS
## A2    A   2            350 80.0000     100          2      0.125             25 1e+02      25       0           640
## B2    B   2            350 40.0000      50          9      0.125             25 1e+02      25      50           320
## C2    C   2            350 20.0000      25         16      0.125             25 1e+02      25      75           160
## D2    D   2            350 10.0000     100         23      0.125             25 1e+02      25       0            80
## E2    E   2            350  5.0000      50         30      0.125             25 1e+02      25      50            40
## F2    F   2            350  2.5000      25         37      0.125             25 1e+02      25      75            20
## G2    G   2            350  1.2500     100         44      0.125             25 1e+02      25       0            10
## H2    H   2            350  0.6250      50         51      0.125             25 1e+02      25      50             5
## I2    I   2            350  0.3125      25         58      0.125             25 1e+02      25      75           2.5
## J2    J   2            350 10.0000     100         65      0.125             25 0e+00       0      25            80
## K2    K   2            350 80.0000     100          2      0.125             25 1e-01      25       0           640
## L2    L   2            350 40.0000      50          9      0.125             25 1e-01      25      50           320
## M2    M   2            350 20.0000      25         16      0.125             25 1e-01      25      75           160
## N2    N   2            350 80.0000     100          2      0.125             25 1e-02      25       0           640
## O2    O   2            350 40.0000      50          9      0.125             25 1e-02      25      50           320
## P2    P   2            350 20.0000      25         16      0.125             25 1e-02      25      75           160
## A3    A   3            350 80.0000     100          3      0.250             25 1e+02      25       0           320
## B3    B   3            350 40.0000      50         10      0.250             25 1e+02      25      50           160
## C3    C   3            350 20.0000      25         17      0.250             25 1e+02      25      75            80
## D3    D   3            350 10.0000     100         24      0.250             25 1e+02      25       0            40
## E3    E   3            350  5.0000      50         31      0.250             25 1e+02      25      50            20
## F3    F   3            350  2.5000      25         38      0.250             25 1e+02      25      75            10
## G3    G   3            350  1.2500     100         45      0.250             25 1e+02      25       0             5
## H3    H   3            350  0.6250      50         52      0.250             25 1e+02      25      50           2.5
## I3    I   3            350  0.3125      25         59      0.250             25 1e+02      25      75          1.25
## J3    J   3            350 10.0000     100         66      0.250             25 0e+00       0      25            40
## K3    K   3            350 80.0000     100          3      0.250             25 1e-01      25       0           320
## L3    L   3            350 40.0000      50         10      0.250             25 1e-01      25      50           160
## M3    M   3            350 20.0000      25         17      0.250             25 1e-01      25      75            80
## N3    N   3            350 80.0000     100          3      0.250             25 1e-02      25       0           320
## O3    O   3            350 40.0000      50         10      0.250             25 1e-02      25      50           160
## P3    P   3            350 20.0000      25         17      0.250             25 1e-02      25      75            80
## A4    A   4            350 80.0000     100          4      0.500             25 1e+02      25       0           160
## B4    B   4            350 40.0000      50         11      0.500             25 1e+02      25      50            80
## C4    C   4            350 20.0000      25         18      0.500             25 1e+02      25      75            40
## D4    D   4            350 10.0000     100         25      0.500             25 1e+02      25       0            20
## E4    E   4            350  5.0000      50         32      0.500             25 1e+02      25      50            10
## F4    F   4            350  2.5000      25         39      0.500             25 1e+02      25      75             5
## G4    G   4            350  1.2500     100         46      0.500             25 1e+02      25       0           2.5
## H4    H   4            350  0.6250      50         53      0.500             25 1e+02      25      50          1.25
## I4    I   4            350  0.3125      25         60      0.500             25 1e+02      25      75         0.625
## J4    J   4            350 10.0000     100         67      0.500             25 0e+00       0      25            20
## K4    K   4            350 80.0000     100          4      0.500             25 1e-01      25       0           160
## L4    L   4            350 40.0000      50         11      0.500             25 1e-01      25      50            80
## M4    M   4            350 20.0000      25         18      0.500             25 1e-01      25      75            40
## N4    N   4            350 80.0000     100          4      0.500             25 1e-02      25       0           160
## O4    O   4            350 40.0000      50         11      0.500             25 1e-02      25      50            80
## P4    P   4            350 20.0000      25         18      0.500             25 1e-02      25      75            40
## A5    A   5            350 80.0000     100          5      1.000             25 1e+02      25       0            80
## B5    B   5            350 40.0000      50         12      1.000             25 1e+02      25      50            40
## C5    C   5            350 20.0000      25         19      1.000             25 1e+02      25      75            20
## D5    D   5            350 10.0000     100         26      1.000             25 1e+02      25       0            10
## E5    E   5            350  5.0000      50         33      1.000             25 1e+02      25      50             5
## F5    F   5            350  2.5000      25         40      1.000             25 1e+02      25      75           2.5
## G5    G   5            350  1.2500     100         47      1.000             25 1e+02      25       0          1.25
## H5    H   5            350  0.6250      50         54      1.000             25 1e+02      25      50         0.625
## I5    I   5            350  0.3125      25         61      1.000             25 1e+02      25      75        0.3125
## J5    J   5            350 10.0000     100         68      1.000             25 0e+00       0      25            10
## K5    K   5            350 80.0000     100          5      1.000             25 1e-01      25       0            80
## L5    L   5            350 40.0000      50         12      1.000             25 1e-01      25      50            40
## M5    M   5            350 20.0000      25         19      1.000             25 1e-01      25      75            20
## N5    N   5            350 80.0000     100          5      1.000             25 1e-02      25       0            80
## O5    O   5            350 40.0000      50         12      1.000             25 1e-02      25      50            40
## P5    P   5            350 20.0000      25         19      1.000             25 1e-02      25      75            20
## A6    A   6            350 80.0000     100          6      2.000             25 1e+02      25       0            40
## B6    B   6            350 40.0000      50         13      2.000             25 1e+02      25      50            20
## C6    C   6            350 20.0000      25         20      2.000             25 1e+02      25      75            10
## D6    D   6            350 10.0000     100         27      2.000             25 1e+02      25       0             5
## E6    E   6            350  5.0000      50         34      2.000             25 1e+02      25      50           2.5
## F6    F   6            350  2.5000      25         41      2.000             25 1e+02      25      75          1.25
## G6    G   6            350  1.2500     100         48      2.000             25 1e+02      25       0         0.625
## H6    H   6            350  0.6250      50         55      2.000             25 1e+02      25      50        0.3125
## I6    I   6            350  0.3125      25         62      2.000             25 1e+02      25      75       0.15625
## J6    J   6            350 10.0000     100         69      2.000             25 0e+00       0      25             5
## K6    K   6            350 80.0000     100          6      2.000             25 1e-01      25       0            40
## L6    L   6            350 40.0000      50         13      2.000             25 1e-01      25      50            20
## M6    M   6            350 20.0000      25         20      2.000             25 1e-01      25      75            10
## N6    N   6            350 80.0000     100          6      2.000             25 1e-02      25       0            40
## O6    O   6            350 40.0000      50         13      2.000             25 1e-02      25      50            20
## P6    P   6            350 20.0000      25         20      2.000             25 1e-02      25      75            10
## A7    A   7            350 80.0000     100          7      4.000             25 1e+02      25       0            20
## B7    B   7            350 40.0000      50         14      4.000             25 1e+02      25      50            10
## C7    C   7            350 20.0000      25         21      4.000             25 1e+02      25      75             5
## D7    D   7            350 10.0000     100         28      4.000             25 1e+02      25       0           2.5
## E7    E   7            350  5.0000      50         35      4.000             25 1e+02      25      50          1.25
## F7    F   7            350  2.5000      25         42      4.000             25 1e+02      25      75         0.625
## G7    G   7            350  1.2500     100         49      4.000             25 1e+02      25       0        0.3125
## H7    H   7            350  0.6250      50         56      4.000             25 1e+02      25      50       0.15625
## I7    I   7            350  0.3125      25         63      4.000             25 1e+02      25      75      0.078125
## J7    J   7            350 10.0000     100         70      4.000             25 0e+00       0      25           2.5
## K7    K   7            350 80.0000     100          7      4.000             25 1e-01      25       0            20
## L7    L   7            350 40.0000      50         14      4.000             25 1e-01      25      50            10
## M7    M   7            350 20.0000      25         21      4.000             25 1e-01      25      75             5
## N7    N   7            350 80.0000     100          7      4.000             25 1e-02      25       0            20
## O7    O   7            350 40.0000      50         14      4.000             25 1e-02      25      50            10
## P7    P   7            350 20.0000      25         21      4.000             25 1e-02      25      75             5
## A8    A   8            350 80.0000     100          1      0.000              0 1e+01      25      25 no_RT_PRIMERS
## B8    B   8            350 40.0000      50          8      0.000              0 1e+01      25      75 no_RT_PRIMERS
## C8    C   8            350 20.0000      25         15      0.000              0 1e+01      25     100 no_RT_PRIMERS
## D8    D   8            350 10.0000     100         22      0.000              0 1e+01      25      25 no_RT_PRIMERS
## E8    E   8            350  5.0000      50         29      0.000              0 1e+01      25      75 no_RT_PRIMERS
## F8    F   8            350  2.5000      25         36      0.000              0 1e+01      25     100 no_RT_PRIMERS
## G8    G   8            350  1.2500     100         43      0.000              0 1e+01      25      25 no_RT_PRIMERS
## H8    H   8            350  0.6250      50         50      0.000              0 1e+01      25      75 no_RT_PRIMERS
## I8    I   8            350  0.3125      25         57      0.000              0 1e+01      25     100 no_RT_PRIMERS
## J8    J   8            350 10.0000     100         64      0.000              0 0e+00       0      50 no_RT_PRIMERS
## K8    K   8            350 10.0000     100         22      0.000              0 1e-01      25      25 no_RT_PRIMERS
## L8    L   8            350  5.0000      50         29      0.000              0 1e-01      25      75 no_RT_PRIMERS
## M8    M   8            350  2.5000      25         36      0.000              0 1e-01      25     100 no_RT_PRIMERS
## N8    N   8            350 10.0000     100         22      0.000              0 1e-02      25      25 no_RT_PRIMERS
## O8    O   8            350  5.0000      50         29      0.000              0 1e-02      25      75 no_RT_PRIMERS
## P8    P   8            350  2.5000      25         36      0.000              0 1e-02      25     100 no_RT_PRIMERS
## A9    A   9            350 80.0000     100          2      0.125             25 1e+01      25       0           640
## B9    B   9            350 40.0000      50          9      0.125             25 1e+01      25      50           320
## C9    C   9            350 20.0000      25         16      0.125             25 1e+01      25      75           160
## D9    D   9            350 10.0000     100         23      0.125             25 1e+01      25       0            80
## E9    E   9            350  5.0000      50         30      0.125             25 1e+01      25      50            40
## F9    F   9            350  2.5000      25         37      0.125             25 1e+01      25      75            20
## G9    G   9            350  1.2500     100         44      0.125             25 1e+01      25       0            10
## H9    H   9            350  0.6250      50         51      0.125             25 1e+01      25      50             5
## I9    I   9            350  0.3125      25         58      0.125             25 1e+01      25      75           2.5
## J9    J   9            350 10.0000     100         65      0.125             25 0e+00       0      25            80
## K9    K   9            350 10.0000     100         23      0.125             25 1e-01      25       0            80
## L9    L   9            350  5.0000      50         30      0.125             25 1e-01      25      50            40
## M9    M   9            350  2.5000      25         37      0.125             25 1e-01      25      75            20
## N9    N   9            350 10.0000     100         23      0.125             25 1e-02      25       0            80
## O9    O   9            350  5.0000      50         30      0.125             25 1e-02      25      50            40
## P9    P   9            350  2.5000      25         37      0.125             25 1e-02      25      75            20
## A10   A  10            350 80.0000     100          3      0.250             25 1e+01      25       0           320
## B10   B  10            350 40.0000      50         10      0.250             25 1e+01      25      50           160
## C10   C  10            350 20.0000      25         17      0.250             25 1e+01      25      75            80
## D10   D  10            350 10.0000     100         24      0.250             25 1e+01      25       0            40
## E10   E  10            350  5.0000      50         31      0.250             25 1e+01      25      50            20
## F10   F  10            350  2.5000      25         38      0.250             25 1e+01      25      75            10
## G10   G  10            350  1.2500     100         45      0.250             25 1e+01      25       0             5
## H10   H  10            350  0.6250      50         52      0.250             25 1e+01      25      50           2.5
## I10   I  10            350  0.3125      25         59      0.250             25 1e+01      25      75          1.25
## J10   J  10            350 10.0000     100         66      0.250             25 0e+00       0      25            40
## K10   K  10            350 10.0000     100         24      0.250             25 1e-01      25       0            40
## L10   L  10            350  5.0000      50         31      0.250             25 1e-01      25      50            20
## M10   M  10            350  2.5000      25         38      0.250             25 1e-01      25      75            10
## N10   N  10            350 10.0000     100         24      0.250             25 1e-02      25       0            40
## O10   O  10            350  5.0000      50         31      0.250             25 1e-02      25      50            20
## P10   P  10            350  2.5000      25         38      0.250             25 1e-02      25      75            10
## A11   A  11            350 80.0000     100          4      0.500             25 1e+01      25       0           160
## B11   B  11            350 40.0000      50         11      0.500             25 1e+01      25      50            80
## C11   C  11            350 20.0000      25         18      0.500             25 1e+01      25      75            40
## D11   D  11            350 10.0000     100         25      0.500             25 1e+01      25       0            20
## E11   E  11            350  5.0000      50         32      0.500             25 1e+01      25      50            10
## F11   F  11            350  2.5000      25         39      0.500             25 1e+01      25      75             5
## G11   G  11            350  1.2500     100         46      0.500             25 1e+01      25       0           2.5
## H11   H  11            350  0.6250      50         53      0.500             25 1e+01      25      50          1.25
## I11   I  11            350  0.3125      25         60      0.500             25 1e+01      25      75         0.625
## J11   J  11            350 10.0000     100         67      0.500             25 0e+00       0      25            20
## K11   K  11            350 10.0000     100         25      0.500             25 1e-01      25       0            20
## L11   L  11            350  5.0000      50         32      0.500             25 1e-01      25      50            10
## M11   M  11            350  2.5000      25         39      0.500             25 1e-01      25      75             5
## N11   N  11            350 10.0000     100         25      0.500             25 1e-02      25       0            20
## O11   O  11            350  5.0000      50         32      0.500             25 1e-02      25      50            10
## P11   P  11            350  2.5000      25         39      0.500             25 1e-02      25      75             5
## A12   A  12            350 80.0000     100          5      1.000             25 1e+01      25       0            80
## B12   B  12            350 40.0000      50         12      1.000             25 1e+01      25      50            40
## C12   C  12            350 20.0000      25         19      1.000             25 1e+01      25      75            20
## D12   D  12            350 10.0000     100         26      1.000             25 1e+01      25       0            10
## E12   E  12            350  5.0000      50         33      1.000             25 1e+01      25      50             5
## F12   F  12            350  2.5000      25         40      1.000             25 1e+01      25      75           2.5
## G12   G  12            350  1.2500     100         47      1.000             25 1e+01      25       0          1.25
## H12   H  12            350  0.6250      50         54      1.000             25 1e+01      25      50         0.625
## I12   I  12            350  0.3125      25         61      1.000             25 1e+01      25      75        0.3125
## J12   J  12            350 10.0000     100         68      1.000             25 0e+00       0      25            10
## K12   K  12            350 10.0000     100         26      1.000             25 1e-01      25       0            10
## L12   L  12            350  5.0000      50         33      1.000             25 1e-01      25      50             5
## M12   M  12            350  2.5000      25         40      1.000             25 1e-01      25      75           2.5
## N12   N  12            350 10.0000     100         26      1.000             25 1e-02      25       0            10
## O12   O  12            350  5.0000      50         33      1.000             25 1e-02      25      50             5
## P12   P  12            350  2.5000      25         40      1.000             25 1e-02      25      75           2.5
## A13   A  13            350 80.0000     100          6      2.000             25 1e+01      25       0            40
## B13   B  13            350 40.0000      50         13      2.000             25 1e+01      25      50            20
## C13   C  13            350 20.0000      25         20      2.000             25 1e+01      25      75            10
## D13   D  13            350 10.0000     100         27      2.000             25 1e+01      25       0             5
## E13   E  13            350  5.0000      50         34      2.000             25 1e+01      25      50           2.5
## F13   F  13            350  2.5000      25         41      2.000             25 1e+01      25      75          1.25
## G13   G  13            350  1.2500     100         48      2.000             25 1e+01      25       0         0.625
## H13   H  13            350  0.6250      50         55      2.000             25 1e+01      25      50        0.3125
## I13   I  13            350  0.3125      25         62      2.000             25 1e+01      25      75       0.15625
## J13   J  13            350 10.0000     100         69      2.000             25 0e+00       0      25             5
## K13   K  13            350 10.0000     100         27      2.000             25 1e-01      25       0             5
## L13   L  13            350  5.0000      50         34      2.000             25 1e-01      25      50           2.5
## M13   M  13            350  2.5000      25         41      2.000             25 1e-01      25      75          1.25
## N13   N  13            350 10.0000     100         27      2.000             25 1e-02      25       0             5
## O13   O  13            350  5.0000      50         34      2.000             25 1e-02      25      50           2.5
## P13   P  13            350  2.5000      25         41      2.000             25 1e-02      25      75          1.25
## A14   A  14            350 80.0000     100          7      4.000             25 1e+01      25       0            20
## B14   B  14            350 40.0000      50         14      4.000             25 1e+01      25      50            10
## C14   C  14            350 20.0000      25         21      4.000             25 1e+01      25      75             5
## D14   D  14            350 10.0000     100         28      4.000             25 1e+01      25       0           2.5
## E14   E  14            350  5.0000      50         35      4.000             25 1e+01      25      50          1.25
## F14   F  14            350  2.5000      25         42      4.000             25 1e+01      25      75         0.625
## G14   G  14            350  1.2500     100         49      4.000             25 1e+01      25       0        0.3125
## H14   H  14            350  0.6250      50         56      4.000             25 1e+01      25      50       0.15625
## I14   I  14            350  0.3125      25         63      4.000             25 1e+01      25      75      0.078125
## J14   J  14            350 10.0000     100         70      4.000             25 0e+00       0      25           2.5
## K14   K  14            350 10.0000     100         28      4.000             25 1e-01      25       0           2.5
## L14   L  14            350  5.0000      50         35      4.000             25 1e-01      25      50          1.25
## M14   M  14            350  2.5000      25         42      4.000             25 1e-01      25      75         0.625
## N14   N  14            350 10.0000     100         28      4.000             25 1e-02      25       0           2.5
## O14   O  14            350  5.0000      50         35      4.000             25 1e-02      25      50          1.25
## P14   P  14            350  2.5000      25         42      4.000             25 1e-02      25      75         0.625
## A15   A  15            350 80.0000     100          1      0.000              0 1e+00      25      25 no_RT_PRIMERS
## B15   B  15            350 40.0000      50          8      0.000              0 1e+00      25      75 no_RT_PRIMERS
## C15   C  15            350 20.0000      25         15      0.000              0 1e+00      25     100 no_RT_PRIMERS
## D15   D  15            350 10.0000     100         22      0.000              0 1e+00      25      25 no_RT_PRIMERS
## E15   E  15            350  5.0000      50         29      0.000              0 1e+00      25      75 no_RT_PRIMERS
## F15   F  15            350  2.5000      25         36      0.000              0 1e+00      25     100 no_RT_PRIMERS
## G15   G  15            350  1.2500     100         43      0.000              0 1e+00      25      25 no_RT_PRIMERS
## H15   H  15            350  0.6250      50         50      0.000              0 1e+00      25      75 no_RT_PRIMERS
## I15   I  15            350  0.3125      25         57      0.000              0 1e+00      25     100 no_RT_PRIMERS
## J15   J  15            350 10.0000     100         64      0.000              0 0e+00       0      50 no_RT_PRIMERS
## K15   K  15            350  1.2500     100         43      0.000              0 1e-01      25      25 no_RT_PRIMERS
## L15   L  15            350  0.6250      50         50      0.000              0 1e-01      25      75 no_RT_PRIMERS
## M15   M  15            350  0.3125      25         57      0.000              0 1e-01      25     100 no_RT_PRIMERS
## N15   N  15            350  1.2500     100         43      0.000              0 1e-02      25      25 no_RT_PRIMERS
## O15   O  15            350  0.6250      50         50      0.000              0 1e-02      25      75 no_RT_PRIMERS
## P15   P  15            350  0.3125      25         57      0.000              0 1e-02      25     100 no_RT_PRIMERS
## A16   A  16            350 80.0000     100          2      0.125             25 1e+00      25       0           640
## B16   B  16            350 40.0000      50          9      0.125             25 1e+00      25      50           320
## C16   C  16            350 20.0000      25         16      0.125             25 1e+00      25      75           160
## D16   D  16            350 10.0000     100         23      0.125             25 1e+00      25       0            80
## E16   E  16            350  5.0000      50         30      0.125             25 1e+00      25      50            40
## F16   F  16            350  2.5000      25         37      0.125             25 1e+00      25      75            20
## G16   G  16            350  1.2500     100         44      0.125             25 1e+00      25       0            10
## H16   H  16            350  0.6250      50         51      0.125             25 1e+00      25      50             5
## I16   I  16            350  0.3125      25         58      0.125             25 1e+00      25      75           2.5
## J16   J  16            350 10.0000     100         65      0.125             25 0e+00       0      25            80
## K16   K  16            350  1.2500     100         44      0.125             25 1e-01      25       0            10
## L16   L  16            350  0.6250      50         51      0.125             25 1e-01      25      50             5
## M16   M  16            350  0.3125      25         58      0.125             25 1e-01      25      75           2.5
## N16   N  16            350  1.2500     100         44      0.125             25 1e-02      25       0            10
## O16   O  16            350  0.6250      50         51      0.125             25 1e-02      25      50             5
## P16   P  16            350  0.3125      25         58      0.125             25 1e-02      25      75           2.5
## A17   A  17            350 80.0000     100          3      0.250             25 1e+00      25       0           320
## B17   B  17            350 40.0000      50         10      0.250             25 1e+00      25      50           160
## C17   C  17            350 20.0000      25         17      0.250             25 1e+00      25      75            80
## D17   D  17            350 10.0000     100         24      0.250             25 1e+00      25       0            40
## E17   E  17            350  5.0000      50         31      0.250             25 1e+00      25      50            20
## F17   F  17            350  2.5000      25         38      0.250             25 1e+00      25      75            10
## G17   G  17            350  1.2500     100         45      0.250             25 1e+00      25       0             5
## H17   H  17            350  0.6250      50         52      0.250             25 1e+00      25      50           2.5
## I17   I  17            350  0.3125      25         59      0.250             25 1e+00      25      75          1.25
## J17   J  17            350 10.0000     100         66      0.250             25 0e+00       0      25            40
## K17   K  17            350  1.2500     100         45      0.250             25 1e-01      25       0             5
## L17   L  17            350  0.6250      50         52      0.250             25 1e-01      25      50           2.5
## M17   M  17            350  0.3125      25         59      0.250             25 1e-01      25      75          1.25
## N17   N  17            350  1.2500     100         45      0.250             25 1e-02      25       0             5
## O17   O  17            350  0.6250      50         52      0.250             25 1e-02      25      50           2.5
## P17   P  17            350  0.3125      25         59      0.250             25 1e-02      25      75          1.25
## A18   A  18            350 80.0000     100          4      0.500             25 1e+00      25       0           160
## B18   B  18            350 40.0000      50         11      0.500             25 1e+00      25      50            80
## C18   C  18            350 20.0000      25         18      0.500             25 1e+00      25      75            40
## D18   D  18            350 10.0000     100         25      0.500             25 1e+00      25       0            20
## E18   E  18            350  5.0000      50         32      0.500             25 1e+00      25      50            10
## F18   F  18            350  2.5000      25         39      0.500             25 1e+00      25      75             5
## G18   G  18            350  1.2500     100         46      0.500             25 1e+00      25       0           2.5
## H18   H  18            350  0.6250      50         53      0.500             25 1e+00      25      50          1.25
## I18   I  18            350  0.3125      25         60      0.500             25 1e+00      25      75         0.625
## J18   J  18            350 10.0000     100         67      0.500             25 0e+00       0      25            20
## K18   K  18            350  1.2500     100         46      0.500             25 1e-01      25       0           2.5
## L18   L  18            350  0.6250      50         53      0.500             25 1e-01      25      50          1.25
## M18   M  18            350  0.3125      25         60      0.500             25 1e-01      25      75         0.625
## N18   N  18            350  1.2500     100         46      0.500             25 1e-02      25       0           2.5
## O18   O  18            350  0.6250      50         53      0.500             25 1e-02      25      50          1.25
## P18   P  18            350  0.3125      25         60      0.500             25 1e-02      25      75         0.625
## A19   A  19            350 80.0000     100          5      1.000             25 1e+00      25       0            80
## B19   B  19            350 40.0000      50         12      1.000             25 1e+00      25      50            40
## C19   C  19            350 20.0000      25         19      1.000             25 1e+00      25      75            20
## D19   D  19            350 10.0000     100         26      1.000             25 1e+00      25       0            10
## E19   E  19            350  5.0000      50         33      1.000             25 1e+00      25      50             5
## F19   F  19            350  2.5000      25         40      1.000             25 1e+00      25      75           2.5
## G19   G  19            350  1.2500     100         47      1.000             25 1e+00      25       0          1.25
## H19   H  19            350  0.6250      50         54      1.000             25 1e+00      25      50         0.625
## I19   I  19            350  0.3125      25         61      1.000             25 1e+00      25      75        0.3125
## J19   J  19            350 10.0000     100         68      1.000             25 0e+00       0      25            10
## K19   K  19            350  1.2500     100         47      1.000             25 1e-01      25       0          1.25
## L19   L  19            350  0.6250      50         54      1.000             25 1e-01      25      50         0.625
## M19   M  19            350  0.3125      25         61      1.000             25 1e-01      25      75        0.3125
## N19   N  19            350  1.2500     100         47      1.000             25 1e-02      25       0          1.25
## O19   O  19            350  0.6250      50         54      1.000             25 1e-02      25      50         0.625
## P19   P  19            350  0.3125      25         61      1.000             25 1e-02      25      75        0.3125
## A20   A  20            350 80.0000     100          6      2.000             25 1e+00      25       0            40
## B20   B  20            350 40.0000      50         13      2.000             25 1e+00      25      50            20
## C20   C  20            350 20.0000      25         20      2.000             25 1e+00      25      75            10
## D20   D  20            350 10.0000     100         27      2.000             25 1e+00      25       0             5
## E20   E  20            350  5.0000      50         34      2.000             25 1e+00      25      50           2.5
## F20   F  20            350  2.5000      25         41      2.000             25 1e+00      25      75          1.25
## G20   G  20            350  1.2500     100         48      2.000             25 1e+00      25       0         0.625
## H20   H  20            350  0.6250      50         55      2.000             25 1e+00      25      50        0.3125
## I20   I  20            350  0.3125      25         62      2.000             25 1e+00      25      75       0.15625
## J20   J  20            350 10.0000     100         69      2.000             25 0e+00       0      25             5
## K20   K  20            350  1.2500     100         48      2.000             25 1e-01      25       0         0.625
## L20   L  20            350  0.6250      50         55      2.000             25 1e-01      25      50        0.3125
## M20   M  20            350  0.3125      25         62      2.000             25 1e-01      25      75       0.15625
## N20   N  20            350  1.2500     100         48      2.000             25 1e-02      25       0         0.625
## O20   O  20            350  0.6250      50         55      2.000             25 1e-02      25      50        0.3125
## P20   P  20            350  0.3125      25         62      2.000             25 1e-02      25      75       0.15625
## A21   A  21            350 80.0000     100          7      4.000             25 1e+00      25       0            20
## B21   B  21            350 40.0000      50         14      4.000             25 1e+00      25      50            10
## C21   C  21            350 20.0000      25         21      4.000             25 1e+00      25      75             5
## D21   D  21            350 10.0000     100         28      4.000             25 1e+00      25       0           2.5
## E21   E  21            350  5.0000      50         35      4.000             25 1e+00      25      50          1.25
## F21   F  21            350  2.5000      25         42      4.000             25 1e+00      25      75         0.625
## G21   G  21            350  1.2500     100         49      4.000             25 1e+00      25       0        0.3125
## H21   H  21            350  0.6250      50         56      4.000             25 1e+00      25      50       0.15625
## I21   I  21            350  0.3125      25         63      4.000             25 1e+00      25      75      0.078125
## J21   J  21            350 10.0000     100         70      4.000             25 0e+00       0      25           2.5
## K21   K  21            350  1.2500     100         49      4.000             25 1e-01      25       0        0.3125
## L21   L  21            350  0.6250      50         56      4.000             25 1e-01      25      50       0.15625
## M21   M  21            350  0.3125      25         63      4.000             25 1e-01      25      75      0.078125
## N21   N  21            350  1.2500     100         49      4.000             25 1e-02      25       0        0.3125
## O21   O  21            350  0.6250      50         56      4.000             25 1e-02      25      50       0.15625
## P21   P  21            350  0.3125      25         63      4.000             25 1e-02      25      75      0.078125
## A22   A  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## B22   B  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## C22   C  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## D22   D  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## E22   E  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## F22   F  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## G22   G  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## H22   H  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## I22   I  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## J22   J  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## K22   K  22            350 10.0000     100         65      0.125             25 0e+00       0      25            80
## L22   L  22            350 10.0000     100         68      1.000             25 0e+00       0      25            10
## M22   M  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## N22   N  22            350 10.0000     100         65      0.125             25 0e+00       0      25            80
## O22   O  22            350 10.0000     100         68      1.000             25 0e+00       0      25            10
## P22   P  22             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## A23   A  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## B23   B  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## C23   C  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## D23   D  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## E23   E  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## F23   F  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## G23   G  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## H23   H  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## I23   I  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## J23   J  23             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## K23   K  23            350 10.0000     100         66      0.250             25 0e+00       0      25            40
## L23   L  23            350 10.0000     100         69      2.000             25 0e+00       0      25             5
## M23   M  23            350 10.0000     100         64      0.000              0 0e+00       0      50 no_RT_PRIMERS
## N23   N  23            350 10.0000     100         66      0.250             25 0e+00       0      25            40
## O23   O  23            350 10.0000     100         69      2.000             25 0e+00       0      25             5
## P23   P  23            350 10.0000     100         64      0.000              0 0e+00       0      50 no_RT_PRIMERS
## A24   A  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## B24   B  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## C24   C  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## D24   D  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## E24   E  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## F24   F  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## G24   G  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## H24   H  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## I24   I  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## J24   J  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## K24   K  24            350 10.0000     100         67      0.500             25 0e+00       0      25            20
## L24   L  24            350 10.0000     100         70      4.000             25 0e+00       0      25           2.5
## M24   M  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
## N24   N  24            350 10.0000     100         67      0.500             25 0e+00       0      25            20
## O24   O  24            350 10.0000     100         70      4.000             25 0e+00       0      25           2.5
## P24   P  24             NA      NA      NA         NA         NA             NA    NA      NA      NA          <NA>
##     total_volume
## A1           500
## B1           500
## C1           500
## D1           500
## E1           500
## F1           500
## G1           500
## H1           500
## I1           500
## J1           500
## K1           500
## L1           500
## M1           500
## N1           500
## O1           500
## P1           500
## A2           500
## B2           500
## C2           500
## D2           500
## E2           500
## F2           500
## G2           500
## H2           500
## I2           500
## J2           500
## K2           500
## L2           500
## M2           500
## N2           500
## O2           500
## P2           500
## A3           500
## B3           500
## C3           500
## D3           500
## E3           500
## F3           500
## G3           500
## H3           500
## I3           500
## J3           500
## K3           500
## L3           500
## M3           500
## N3           500
## O3           500
## P3           500
## A4           500
## B4           500
## C4           500
## D4           500
## E4           500
## F4           500
## G4           500
## H4           500
## I4           500
## J4           500
## K4           500
## L4           500
## M4           500
## N4           500
## O4           500
## P4           500
## A5           500
## B5           500
## C5           500
## D5           500
## E5           500
## F5           500
## G5           500
## H5           500
## I5           500
## J5           500
## K5           500
## L5           500
## M5           500
## N5           500
## O5           500
## P5           500
## A6           500
## B6           500
## C6           500
## D6           500
## E6           500
## F6           500
## G6           500
## H6           500
## I6           500
## J6           500
## K6           500
## L6           500
## M6           500
## N6           500
## O6           500
## P6           500
## A7           500
## B7           500
## C7           500
## D7           500
## E7           500
## F7           500
## G7           500
## H7           500
## I7           500
## J7           500
## K7           500
## L7           500
## M7           500
## N7           500
## O7           500
## P7           500
## A8           500
## B8           500
## C8           500
## D8           500
## E8           500
## F8           500
## G8           500
## H8           500
## I8           500
## J8           500
## K8           500
## L8           500
## M8           500
## N8           500
## O8           500
## P8           500
## A9           500
## B9           500
## C9           500
## D9           500
## E9           500
## F9           500
## G9           500
## H9           500
## I9           500
## J9           500
## K9           500
## L9           500
## M9           500
## N9           500
## O9           500
## P9           500
## A10          500
## B10          500
## C10          500
## D10          500
## E10          500
## F10          500
## G10          500
## H10          500
## I10          500
## J10          500
## K10          500
## L10          500
## M10          500
## N10          500
## O10          500
## P10          500
## A11          500
## B11          500
## C11          500
## D11          500
## E11          500
## F11          500
## G11          500
## H11          500
## I11          500
## J11          500
## K11          500
## L11          500
## M11          500
## N11          500
## O11          500
## P11          500
## A12          500
## B12          500
## C12          500
## D12          500
## E12          500
## F12          500
## G12          500
## H12          500
## I12          500
## J12          500
## K12          500
## L12          500
## M12          500
## N12          500
## O12          500
## P12          500
## A13          500
## B13          500
## C13          500
## D13          500
## E13          500
## F13          500
## G13          500
## H13          500
## I13          500
## J13          500
## K13          500
## L13          500
## M13          500
## N13          500
## O13          500
## P13          500
## A14          500
## B14          500
## C14          500
## D14          500
## E14          500
## F14          500
## G14          500
## H14          500
## I14          500
## J14          500
## K14          500
## L14          500
## M14          500
## N14          500
## O14          500
## P14          500
## A15          500
## B15          500
## C15          500
## D15          500
## E15          500
## F15          500
## G15          500
## H15          500
## I15          500
## J15          500
## K15          500
## L15          500
## M15          500
## N15          500
## O15          500
## P15          500
## A16          500
## B16          500
## C16          500
## D16          500
## E16          500
## F16          500
## G16          500
## H16          500
## I16          500
## J16          500
## K16          500
## L16          500
## M16          500
## N16          500
## O16          500
## P16          500
## A17          500
## B17          500
## C17          500
## D17          500
## E17          500
## F17          500
## G17          500
## H17          500
## I17          500
## J17          500
## K17          500
## L17          500
## M17          500
## N17          500
## O17          500
## P17          500
## A18          500
## B18          500
## C18          500
## D18          500
## E18          500
## F18          500
## G18          500
## H18          500
## I18          500
## J18          500
## K18          500
## L18          500
## M18          500
## N18          500
## O18          500
## P18          500
## A19          500
## B19          500
## C19          500
## D19          500
## E19          500
## F19          500
## G19          500
## H19          500
## I19          500
## J19          500
## K19          500
## L19          500
## M19          500
## N19          500
## O19          500
## P19          500
## A20          500
## B20          500
## C20          500
## D20          500
## E20          500
## F20          500
## G20          500
## H20          500
## I20          500
## J20          500
## K20          500
## L20          500
## M20          500
## N20          500
## O20          500
## P20          500
## A21          500
## B21          500
## C21          500
## D21          500
## E21          500
## F21          500
## G21          500
## H21          500
## I21          500
## J21          500
## K21          500
## L21          500
## M21          500
## N21          500
## O21          500
## P21          500
## A22           NA
## B22           NA
## C22           NA
## D22           NA
## E22           NA
## F22           NA
## G22           NA
## H22           NA
## I22           NA
## J22           NA
## K22          500
## L22          500
## M22           NA
## N22          500
## O22          500
## P22           NA
## A23           NA
## B23           NA
## C23           NA
## D23           NA
## E23           NA
## F23           NA
## G23           NA
## H23           NA
## I23           NA
## J23           NA
## K23          500
## L23          500
## M23          500
## N23          500
## O23          500
## P23          500
## A24           NA
## B24           NA
## C24           NA
## D24           NA
## E24           NA
## F24           NA
## G24           NA
## H24           NA
## I24           NA
## J24           NA
## K24          500
## L24          500
## M24           NA
## N24          500
## O24          500
## P24           NA
```

```r
plate$total_volume
```

```
##   [1] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
##  [29] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
##  [57] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
##  [85] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [113] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [141] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [169] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [197] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [225] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [253] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [281] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [309] 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500
## [337]  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA 500 500  NA 500 500  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA 500 500
## [365] 500 500 500 500  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA 500 500  NA 500 500  NA
```

```r
plate$BARCODE_ID
```

```
##   [1]  1  8 15 22 29 36 43 50 57 64  1  8 15  1  8 15  2  9 16 23 30 37 44 51 58 65  2  9 16  2  9 16  3 10 17 24 31 38
##  [39] 45 52 59 66  3 10 17  3 10 17  4 11 18 25 32 39 46 53 60 67  4 11 18  4 11 18  5 12 19 26 33 40 47 54 61 68  5 12
##  [77] 19  5 12 19  6 13 20 27 34 41 48 55 62 69  6 13 20  6 13 20  7 14 21 28 35 42 49 56 63 70  7 14 21  7 14 21  1  8
## [115] 15 22 29 36 43 50 57 64 22 29 36 22 29 36  2  9 16 23 30 37 44 51 58 65 23 30 37 23 30 37  3 10 17 24 31 38 45 52
## [153] 59 66 24 31 38 24 31 38  4 11 18 25 32 39 46 53 60 67 25 32 39 25 32 39  5 12 19 26 33 40 47 54 61 68 26 33 40 26
## [191] 33 40  6 13 20 27 34 41 48 55 62 69 27 34 41 27 34 41  7 14 21 28 35 42 49 56 63 70 28 35 42 28 35 42  1  8 15 22
## [229] 29 36 43 50 57 64 43 50 57 43 50 57  2  9 16 23 30 37 44 51 58 65 44 51 58 44 51 58  3 10 17 24 31 38 45 52 59 66
## [267] 45 52 59 45 52 59  4 11 18 25 32 39 46 53 60 67 46 53 60 46 53 60  5 12 19 26 33 40 47 54 61 68 47 54 61 47 54 61
## [305]  6 13 20 27 34 41 48 55 62 69 48 55 62 48 55 62  7 14 21 28 35 42 49 56 63 70 49 56 63 49 56 63 NA NA NA NA NA NA
## [343] NA NA NA NA 65 68 NA 65 68 NA NA NA NA NA NA NA NA NA NA NA 66 69 64 66 69 64 NA NA NA NA NA NA NA NA NA NA 67 70
## [381] NA 67 70 NA
```

```r
length(which(plate$total_volume != "NA")) 
```

```
## [1] 350
```

```r
sum(is.na(plate$total_volume))
```

```
## [1] 34
```

```r
count(plate$BARCODE_ID)
```

```
##     x freq
## 1   1    5
## 2   2    5
## 3   3    5
## 4   4    5
## 5   5    5
## 6   6    5
## 7   7    5
## 8   8    5
## 9   9    5
## 10 10    5
## 11 11    5
## 12 12    5
## 13 13    5
## 14 14    5
## 15 15    5
## 16 16    5
## 17 17    5
## 18 18    5
## 19 19    5
## 20 20    5
## 21 21    5
## 22 22    5
## 23 23    5
## 24 24    5
## 25 25    5
## 26 26    5
## 27 27    5
## 28 28    5
## 29 29    5
## 30 30    5
## 31 31    5
## 32 32    5
## 33 33    5
## 34 34    5
## 35 35    5
## 36 36    5
## 37 37    5
## 38 38    5
## 39 39    5
## 40 40    5
## 41 41    5
## 42 42    5
## 43 43    5
## 44 44    5
## 45 45    5
## 46 46    5
## 47 47    5
## 48 48    5
## 49 49    5
## 50 50    5
## 51 51    5
## 52 52    5
## 53 53    5
## 54 54    5
## 55 55    5
## 56 56    5
## 57 57    5
## 58 58    5
## 59 59    5
## 60 60    5
## 61 61    5
## 62 62    5
## 63 63    5
## 64 64    5
## 65 65    5
## 66 66    5
## 67 67    5
## 68 68    5
## 69 69    5
## 70 70    5
## 71 NA   34
```

```r
plate$PRIMERS_RATIO
```

```
##   [1] "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS"
##   [8] "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS"
##  [15] "no_RT_PRIMERS" "no_RT_PRIMERS" "640"           "320"           "160"           "80"            "40"           
##  [22] "20"            "10"            "5"             "2.5"           "80"            "640"           "320"          
##  [29] "160"           "640"           "320"           "160"           "320"           "160"           "80"           
##  [36] "40"            "20"            "10"            "5"             "2.5"           "1.25"          "40"           
##  [43] "320"           "160"           "80"            "320"           "160"           "80"            "160"          
##  [50] "80"            "40"            "20"            "10"            "5"             "2.5"           "1.25"         
##  [57] "0.625"         "20"            "160"           "80"            "40"            "160"           "80"           
##  [64] "40"            "80"            "40"            "20"            "10"            "5"             "2.5"          
##  [71] "1.25"          "0.625"         "0.3125"        "10"            "80"            "40"            "20"           
##  [78] "80"            "40"            "20"            "40"            "20"            "10"            "5"            
##  [85] "2.5"           "1.25"          "0.625"         "0.3125"        "0.15625"       "5"             "40"           
##  [92] "20"            "10"            "40"            "20"            "10"            "20"            "10"           
##  [99] "5"             "2.5"           "1.25"          "0.625"         "0.3125"        "0.15625"       "0.078125"     
## [106] "2.5"           "20"            "10"            "5"             "20"            "10"            "5"            
## [113] "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS"
## [120] "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS"
## [127] "no_RT_PRIMERS" "no_RT_PRIMERS" "640"           "320"           "160"           "80"            "40"           
## [134] "20"            "10"            "5"             "2.5"           "80"            "80"            "40"           
## [141] "20"            "80"            "40"            "20"            "320"           "160"           "80"           
## [148] "40"            "20"            "10"            "5"             "2.5"           "1.25"          "40"           
## [155] "40"            "20"            "10"            "40"            "20"            "10"            "160"          
## [162] "80"            "40"            "20"            "10"            "5"             "2.5"           "1.25"         
## [169] "0.625"         "20"            "20"            "10"            "5"             "20"            "10"           
## [176] "5"             "80"            "40"            "20"            "10"            "5"             "2.5"          
## [183] "1.25"          "0.625"         "0.3125"        "10"            "10"            "5"             "2.5"          
## [190] "10"            "5"             "2.5"           "40"            "20"            "10"            "5"            
## [197] "2.5"           "1.25"          "0.625"         "0.3125"        "0.15625"       "5"             "5"            
## [204] "2.5"           "1.25"          "5"             "2.5"           "1.25"          "20"            "10"           
## [211] "5"             "2.5"           "1.25"          "0.625"         "0.3125"        "0.15625"       "0.078125"     
## [218] "2.5"           "2.5"           "1.25"          "0.625"         "2.5"           "1.25"          "0.625"        
## [225] "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS"
## [232] "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS" "no_RT_PRIMERS"
## [239] "no_RT_PRIMERS" "no_RT_PRIMERS" "640"           "320"           "160"           "80"            "40"           
## [246] "20"            "10"            "5"             "2.5"           "80"            "10"            "5"            
## [253] "2.5"           "10"            "5"             "2.5"           "320"           "160"           "80"           
## [260] "40"            "20"            "10"            "5"             "2.5"           "1.25"          "40"           
## [267] "5"             "2.5"           "1.25"          "5"             "2.5"           "1.25"          "160"          
## [274] "80"            "40"            "20"            "10"            "5"             "2.5"           "1.25"         
## [281] "0.625"         "20"            "2.5"           "1.25"          "0.625"         "2.5"           "1.25"         
## [288] "0.625"         "80"            "40"            "20"            "10"            "5"             "2.5"          
## [295] "1.25"          "0.625"         "0.3125"        "10"            "1.25"          "0.625"         "0.3125"       
## [302] "1.25"          "0.625"         "0.3125"        "40"            "20"            "10"            "5"            
## [309] "2.5"           "1.25"          "0.625"         "0.3125"        "0.15625"       "5"             "0.625"        
## [316] "0.3125"        "0.15625"       "0.625"         "0.3125"        "0.15625"       "20"            "10"           
## [323] "5"             "2.5"           "1.25"          "0.625"         "0.3125"        "0.15625"       "0.078125"     
## [330] "2.5"           "0.3125"        "0.15625"       "0.078125"      "0.3125"        "0.15625"       "0.078125"     
## [337] NA              NA              NA              NA              NA              NA              NA             
## [344] NA              NA              NA              "80"            "10"            NA              "80"           
## [351] "10"            NA              NA              NA              NA              NA              NA             
## [358] NA              NA              NA              NA              NA              "40"            "5"            
## [365] "no_RT_PRIMERS" "40"            "5"             "no_RT_PRIMERS" NA              NA              NA             
## [372] NA              NA              NA              NA              NA              NA              NA             
## [379] "20"            "2.5"           NA              "20"            "2.5"           NA
```
