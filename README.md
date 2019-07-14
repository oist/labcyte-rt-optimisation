Optimisation of reverse-transcription using the Labcyte Echo platform
=====================================================================

With a Labcyte Echo 525 instrument, we generated a large number of combinations
of reaction parameters for reverse-transcription in a 384-well plate.  We
read the results by preparing nanoCAGE libraries of the cDNAs.

See also <http://dgt-gitlab.gsc.riken.jp/gitlab/plessy/Labcyte-qPCR> for a proof
of principle on qPCR optimisation.


Transfer designs
----------------

 - [1](Labcyte-RT.md): [TSO] vs [RT] vs [RNA]
 - [2](Labcyte-RT2.md): similar to 2 with different concentration ranges.
 - [3](Labcyte-RT3.md): similar to 3 with different wells used for some TSOs.
 - [4](Labcyte-RT5.md) and [5](Labcyte-RT5.md): fixed concentrations, randomised positions for the TSOs.
 - [6a](Labcyte-RT6a.md), [6b](Labcyte-RT6b.md), [6c](Labcyte-RT6c.md) and [6d](Labcyte-RT6d.md):
   [TSO] vs [RT] vs [RNA] with randomised positions.
 - [7](Labcyte-RT7.md): [dNTP] vs [Mg<sup>2+</sup>] vs [Mn<sup>2+</sup>] with randomised positions.


Concentration checks
--------------------

 - [1](TSO_concentration_check.md): source plates 1 and 2.  Confirmed that the TSOs
   in "row B" were present at the correct concentration.
 - [2](TSO_concentration_check2.md): "old" stock plate (PO_8268526).
   Showed variations of concentration.
 - [3](TSO_concentration_check3.md): source plate for Experiment 6.
   Suggested that attempt to correct variations of concentration in the "old"
   stock plate did not have effect.
 - [4](TSO_concentration_check4.md): source plate for Experiment 7.  Stock TSOs
   were diluted assuming a molarity of 1 mM.  Concentration check number
   [2](TSO_concentration_check2.md) showed that the real molarity was a bit lower.
   The measurements here show that this was reflected in the source plate,
   except for the highest concentrations.
 - [5](TSO_concentration_check5.md): source plate for Experiment 8.  Confirms that
   concentration corrections were effective, but that there is a systematic
   shift (~15 µM instead of 10 as planned).

Data
----

 - nanoCAGE sequence data: [Zenodo 1680999](https://doi.org/10.5281/zenodo.1680999)
 - nanoCAGE sequence alignments: [Zenodo 1683162](https://doi.org/10.5281/zenodo.1683162)


Experiments
-----------
 
### Experiment 1 (pilot) and 2 (replicate)
 
[Experiment one](Labcyte-RT_Data_Analysis.md): plate `A`.
MiSeq ID: `171227_M00528_0321_000000000-B4GLP` (moderate quality, especially on read 1).
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/171227_M00528_0321_000000000-B4GLP.paired_raw_quality_control2.20171228143720/171227_M00528_0321_000000000-B4GLP.paired_raw_quality_control2.20171228143720.html),
  [Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/171227_M00528_0321_000000000-B4GLP.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180104125850/171227_M00528_0321_000000000-B4GLP.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180104125850.html).

Transfer plate design [1](Labcyte-RT.md) was used to compare [TSO] vs [RT] vs [RNA]
in a single 384-well plate.  RTs made on the same RNA amounts were multiplexed
in libraries amplified separately.  Only one amplification was successful.  Analysis
suggested that reaction efficiency reaches a maximum at 20 μM TSO.  It also
showed an obvious artefact at 40 μM TSO, where reaction efficiency was abnormally low.

Barcode extraction was assessed by comparing TagDust 2 with a simple match of
barcode sequence (no HMMs nor correction).  For most of them, the concordance
was very high.  For a few of them, TagDust 2 detected more, probably thanks
to error correction.  These barcodes were not related to the artefact at 40 μM TSO.

Explanation for the artefact was found later, but this experiment was not re-analysed
(superseded by Experiemnt 2 anyway).


[Experiment two](Labcyte-RT_Data_Analysis_2.md): plates: `B`, `C`, and `D`.
MiSeq ID: `180123_M00528_0325_000000000-B4PCK` (bad quality, especially on read 1, leading to large numbers of unmapped reads). 
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180123_M00528_0325_000000000-B4PCK.paired_raw_quality_control2.20180124101336/180123_M00528_0325_000000000-B4PCK.paired_raw_quality_control2.20180124101336.html),
  [Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180123_M00528_0325_000000000-B4PCK.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180124102551/180123_M00528_0325_000000000-B4PCK.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180124102551.html).

Same transfer design ([1](Labcyte-RT.md)) was used to repeat Experiment 1,
this time in triplicates.  All observations, including the artefact, were
confirmed.

Later re-analysis with the correction factors calculated in Experiment 5 suggested
that highest TSO concentrations give the highest yields and lowest amounts
of rDNA and high promoter rates.


### Experiment 3 and 4 (different layouts of the source plate)

[Experiment three](Labcyte-RT_Data_Analysis_3merge.md)
   ([SSIII QC](Labcyte-RT_Data_Analysis_3a.md): plates `E`, `F` and `G`,
   [SSIV QC](Labcyte-RT_Data_Analysis_3b.md): plates `H`, `I` and `J`).
MiSeq ID: `180403_M00528_0348_000000000-B4GP8`.
MOIRAI
 [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180403_M00528_0348_000000000-B4GP8.paired_raw_quality_control2.20180404155647/180403_M00528_0348_000000000-B4GP8.paired_raw_quality_control2.20180404155647.html),
 [Workflow plates 1~3](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180403_M00528_0348_000000000-B4GP8_p123.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180409105009/180403_M00528_0348_000000000-B4GP8_p123.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180409105009.html),
 [Workflow plates 4~6](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180403_M00528_0348_000000000-B4GP8_p456.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180411094125/180403_M00528_0348_000000000-B4GP8_p456.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180411094125.html).

A new transfer design ([2](Labcyte-RT2.md)) was used to repeat Experiment 2,
in triplicates again, this time with two different enzymes.  The artefact was
shown to be related to TSOs in source plate row B, and not related to the
actual concentration.


[Experiment four](Labcyte-RT_Data_Analysis_4.md): plates `K` and `L`.
MiSeq ID: `180326_M00528_0346_000000000-B4GJR`.
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180326_M00528_0346_000000000-B4GJR.paired_raw_quality_control2.20180403111017/180326_M00528_0346_000000000-B4GJR.paired_raw_quality_control2.20180403111017.html),
  [Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180326_M00528_0346_000000000-B4GJR.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180329132046/180326_M00528_0346_000000000-B4GJR.OP-WORKFLOW-CAGEscan-short-reads-v2.1~rc1.20180329132046.html).
     
A new transfer design ([3](Labcyte-RT3.md)) was used to change of position
or contents of "row B" TSOs.  The conclusion was that the artefact was related
to the TSO quality or barcode sequence.


### Experiment 5 and 6 (showing batch bias of the TSOs)
     
[Experiment five](Labcyte-RT_Data_Analysis_5.md): plates `M` and `N`.
MiSeq ID: `180411_M00528_0351_000000000-BN3BL`.
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180411_M00528_0351_000000000-BN3BL.paired_raw_quality_control2.20180412055341/180411_M00528_0351_000000000-BN3BL.paired_raw_quality_control2.20180412055341.html),
  [Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180411_M00528_0351_000000000-BN3BL.OP-WORKFLOW-CAGEscan-short-reads-v2.1.1.20180412203518/180411_M00528_0351_000000000-BN3BL.OP-WORKFLOW-CAGEscan-short-reads-v2.1.1.20180412203518.html).

A new transfer design with randomised TSO positions ([4](Labcyte-RT4.md) and
[5](Labcyte-RT5.md) for two independant randomisations) showed a dramatic
variation of reaction efficiencies, even when the only degree of freedom is the
TSO.  Thus, the origin of the artefact is either the barcode sequence or the
the batch of synthesis.

From these results, normalisation factors were calculated, to correct for
TSO-specific variations of reaction efficiency.


[Experiment six](Labcyte-RT_Data_Analysis_6.md): plates `O` and `P`.
MiSeq ID: `180501_M00528_0359_000000000-B4PJY`.
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180501_M00528_0359_000000000-B4PJY.paired_raw_quality_control2.20180502080443/180501_M00528_0359_000000000-B4PJY.paired_raw_quality_control2.20180502080443.html),
  [Workflow](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180501_M00528_0359_000000000-B4PJY.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180502081812/180501_M00528_0359_000000000-B4PJY.OP-WORKFLOW-CAGEscan-short-reads-v2.1.2.20180502081812.html).

Same randomised transfer designs ([4](Labcyte-RT4.md) and
[5](Labcyte-RT5.md)), but a different batch of TSOs.  In this experiment, the
variation of reaction efficiency is much milder, demonstrating that the problem
was the quality of the TSOs used in the previous experiments.


### Experiment 7 (SSIII) and 8 (SSIV)

Experiment seven
  ([Analysis](Labcyte-RT_Data_Analysis_7.md)) /
  ([QC](Labcyte-RT_Data_Analysis_7_QC.md)): plates `Q`, `R`, `S` and `T`.
MiSeq ID `180517_M00528_0364_000000000-BRGK6` (good quality).
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180517_M00528_0364_000000000-BRGK6.paired_raw_quality_control2.20180518043259/180517_M00528_0364_000000000-BRGK6.paired_raw_quality_control2.20180518043259.html),
  [Workflow]().

Experiment eight, using the same design as exp. 7 but with SuperScript IV.
  ([Analysis](Labcyte-RT_Data_Analysis_8.md)) /
  ([QC](Labcyte-RT_Data_Analysis_8_QC.md)): plates `Q2`, `R2`, `S2` and `T2`.
MiSeq ID `180606_M00528_0367_000000000-BN3FG` (good quality).
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180606_M00528_0367_000000000-BN3FG.paired_raw_quality_control2.20180607095712/180606_M00528_0367_000000000-BN3FG.paired_raw_quality_control2.20180607095712.html),
  [Workflow]().


### Experiment 9 (failed)

Experiment nine, testing dNTP and divalent ion concentration using [transfer
design 7](Labcyte-RT7.md).
  ([Analysis](Labcyte-RT_Data_Analysis_9.md)) /
  ([QC](Labcyte-RT_Data_Analysis_9_QC.md)): plates `U`, `V`, `W` and `X`.
MiSeq ID `180607_M00528_0368_000000000-BN9KM` (good sequence quality, but few reads).
MOIRAI
  [QC](http://moirai.gsc.riken.jp/osc-fs_home/scratch/moirai/nanoCAGE2/project/Labcyte/180607_M00528_0368_000000000-BN9KM.paired_raw_quality_control2.20180608154048/180607_M00528_0368_000000000-BN9KM.paired_raw_quality_control2.20180608154048.html),
  [Workflow]().


### Single cells

Our results match well with a collection of [experiments on single
cells](singleCells.md) that we performed earlier.  The raw data is available on
Zenodo (<http://doi.org/10.5281/zenodo.250156>).  The alignments will be
uploaded soon.


Xmas project
------------

 - [Xmas project](Xmas.md)
