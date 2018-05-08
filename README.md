Optimisation of reverse-transcription using the Labcyte Echo platform
=====================================================================

With the Labcyte Echo platforms, we generate a large number of combinations
of reaction parameters for reverse-transcription in a 384-well plate.  We
read the results by preparing nanoCAGE libraries of the cDNAs.

See also <http://dgt-gitlab.gsc.riken.jp/gitlab/plessy/Labcyte-qPCR> for a proof
of principle on qPCR optimisation.

Source plates designs
---------------------

 - [1](Labcyte-RT.md): [TSO] vs [RT] vs [RNA]
 - [2](Labcyte-RT2.md): similar to 2 with different concentration ranges.
 - [3](Labcyte-RT3.md): similar to 3 with different wells used for some TSOs.
 - [4](Labcyte-RT5.md) and [5](Labcyte-RT5.md): fixed concentrations, randomised positions for the TSOs.
 - [6](Labcyte-RT6.md]): [TSO] vs [RT] vs [RNA] with randomised positions.

 - [Concentration checks](TSO_concentration_check.md) in source plates 1 and 2.

Experiments
-----------
 
### Experiment 1
 
[Experiment one](Labcyte-RT_Data_Analysis.md)
MiSeq ID: `171227_M00528_0321_000000000-B4GLP`.
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/171227_M00528_0321_000000000-B4GLP.paired_raw_quality_control2.20171228143720/171227_M00528_0321_000000000-B4GLP.paired_raw_quality_control2.20171228143720.html),
  [Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/171227_M00528_0321_000000000-B4GLP.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180104125850/171227_M00528_0321_000000000-B4GLP.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180104125850.html).

Source plate design [1](Labcyte-RT.md) was used to compare [TSO] vs [RT] vs [RNA]
in a single 384-well plate.  RTs made on the same RNA amounts were multiplexed
in libraries amplified separately.  Only one amplification was successful.  Analysis
suggested that reaction efficiency reaches a maximum at 20 μM TSO.  It also
showed an obvious artefact at 20 μM TSO, where reaction efficiency was abnormally low.

Explanation for the artefact was found later, but this experiment was not re-analysed
(superseded by Experiemnt 2 anyway).

### Experiment 2

[Experiment two](Labcyte-RT_Data_Analysis_2.md) (same model, three replicates). MiSeq ID: `180123_M00528_0325_000000000-B4PCK`. MOIRAI
     [[QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180123_M00528_0325_000000000-B4PCK.paired_raw_quality_control2.20180124101336/180123_M00528_0325_000000000-B4PCK.paired_raw_quality_control2.20180124101336.html)]
     [[Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180123_M00528_0325_000000000-B4PCK.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180124102551/180123_M00528_0325_000000000-B4PCK.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180124102551.html)]

### Experiment 3

[Experiment three](Labcyte-RT_Data_Analysis_3merge.md)
   [[SSIII QC](Labcyte-RT_Data_Analysis_3a.md)]
   [[SSIV QC](Labcyte-RT_Data_Analysis_3b.md)]]
   (same plate model, different concentrations, three replicates, two enzymes).
   MiSeq ID: `180403_M00528_0348_000000000-B4GP8`.
   MOIRAI [[QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180403_M00528_0348_000000000-B4GP8.paired_raw_quality_control2.20180404155647/180403_M00528_0348_000000000-B4GP8.paired_raw_quality_control2.20180404155647.html)]
     [[Workflow plates 1~3]](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180403_M00528_0348_000000000-B4GP8_p123.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180409105009/180403_M00528_0348_000000000-B4GP8_p123.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180409105009.html)
     [[Workflow plates 4~6]](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180403_M00528_0348_000000000-B4GP8_p456.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180411094125/180403_M00528_0348_000000000-B4GP8_p456.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180411094125.html)

### Experiment 4

[Experiment four](Labcyte-RT_Data_Analysis_4.md) (change of position or contents of "row B" TSOs). MiSeq ID: `180326_M00528_0346_000000000-B4GJR`. MOIRAI
     [[QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180326_M00528_0346_000000000-B4GJR.paired_raw_quality_control2.20180403111017/180326_M00528_0346_000000000-B4GJR.paired_raw_quality_control2.20180403111017.html)]
     [[Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180326_M00528_0346_000000000-B4GJR.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180329132046/180326_M00528_0346_000000000-B4GJR.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180329132046.html)]
     
### Experiment 5
     
     [Experiment five](Labcyte-RT_Data_Analysis_5.md) (Randomisation of TSO position). MiSeq ID: `180411_M00528_0351_000000000-BN3BL`. MOIRAI
     [[QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180411_M00528_0351_000000000-BN3BL.paired_raw_quality_control2.20180412055341/180411_M00528_0351_000000000-BN3BL.paired_raw_quality_control2.20180412055341.html)]
     [[Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180411_M00528_0351_000000000-BN3BL.OP-WORKFLOW-CAGEscan-short-reads-v2.1.1.20180412203518/180411_M00528_0351_000000000-BN3BL.OP-WORKFLOW-CAGEscan-short-reads-v2.1.1.20180412203518.html)]

### Experiment 6

[Experiment six](Labcyte-RT_Data_Analysis_6.md) (Randomisation of TSO position, older TSO master plate). MiSeq ID: `180501_M00528_0359_000000000-B4PJY`. MOIRAI
     [[QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180501_M00528_0359_000000000-B4PJY.paired_raw_quality_control2.20180502080443/180501_M00528_0359_000000000-B4PJY.paired_raw_quality_control2.20180502080443.html)]
     [[Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180501_M00528_0359_000000000-B4PJY.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180502081812/180501_M00528_0359_000000000-B4PJY.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180502081812.html)]

Xmas project
------------

 - [Xmas project](Xmas.md)
