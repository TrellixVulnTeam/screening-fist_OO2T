# Screening Fist

<div class="md-source__icon md-icon">
<a href="https://github.com/jamesengleback/screening-fist">
<svg viewBox="0 0 448 512" xmlns="http://www.w3.org/2000/svg"><!--! Font Awesome Free 6.1.1 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license/free (Icons: CC BY 4.0, Fonts: SIL OFL 1.1, Code: MIT License) Copyright 2022 Fonticons, Inc.--><path d="M439.55 236.05 244 40.45a28.87 28.87 0 0 0-40.81 0l-40.66 40.63 51.52 51.52c27.06-9.14 52.68 16.77 43.39 43.68l49.66 49.66c34.23-11.8 61.18 31 35.47 56.69-26.49 26.49-70.21-2.87-56-37.34L240.22 199v121.85c25.3 12.54 22.26 41.85 9.08 55a34.34 34.34 0 0 1-48.55 0c-17.57-17.6-11.07-46.91 11.25-56v-123c-20.8-8.51-24.6-30.74-18.64-45L142.57 101 8.45 235.14a28.86 28.86 0 0 0 0 40.81l195.61 195.6a28.86 28.86 0 0 0 40.8 0l194.69-194.69a28.86 28.86 0 0 0 0-40.81z"></path></svg>
</a>
</div> 
### [Github](https://github.com/jamesengleback/screening-fist)

## About

Screening Fist is an enzyme:compound screening program focussed on mutants of the P450 BM3 and pharamceutical compounds.
Binding data generated here is used to train a predictive binding model which itself is used to design subsequent screening rounds based on expected information gain.

## Contents

- [Assay Development](assay-development.md)
- [Enzymes Screening Panel](enzymes.md)
- [Compounds Screening Panel](compounds.md)
- [Protocol](protocol.md)
- [Screen Design](screen-design.md)
- [Data Processing and Analysis](data.md)
- [Data Mining](data-mining.md)
- [Model and Training](nn.md)
- [Model-based Experimental Design](al.md)

## Todo

```mermaid
gantt
        title Schedule
        %dateFormat YYYY-MM-DD
        %%axisFormat %Y/%m/%d
        %%axisFormat  %Y-%m-%d


        Deadline                                : milestone deadline, 2022-06-30, 1d

        section Model
        Train Prototype                      : sxfst0, 3d
        Dropout Inference                    : sxfst1, 3d
	Re-Train on Screening data 		: sxfst3, after sxfst12, 2d
	Evaluate 			: sxfst4, after sxfst3, 2d
	Design Next Experiment 		: sxfst5, after sxfst4, 2d
	Methods 				: sxfst2, 14d

	section Screening Data
	Michaelis Menten Metrics 	     : sxfst10, after sxfst1, 2d
	Hit/Miss detection 	     : sxfst11, after sxfst1, 2d
	Make Dataset    	     : sxfst12, after sxfst11, 2d
	Methods 			:  sxfst12, 7d

	section Results
	Pre-Training Data Analysis 			:  sxfst20, 7d
	Screening Data Analysis 			:  sxfst21, after sxfst12, 7d
	Model Analysis 			: sxfst21, after sxfst4, 7d
	Experiment Design Analysis 	: sxfst22, after sxfst5, 7d



```

## Commits

``` git
*   commit ca7fa815d9cd948179b9631a4fa6bb187a02b163
|\  Merge: a3f124c 723f9b9
| | Author: jamesengleback <jamesengleback@hotmail.co.uk>
| | Date:   Sat Jun 4 09:10:38 2022 +0100
| | 
| |     Merge branch 'master' of https://github.com/jamesengleback/screening-fist
| |     merging model changes with docs changes
| | 
| * commit 723f9b97cf6e3413b4b80f2d91b6f2b9928a50cf
| | Author: jamesengleback <jamesengleback@hotmail.co.uk>
| | Date:   Sat Jun 4 07:21:31 2022 +0000
| | 
| |     train.py: fewer saves - disk space. syntax error fix
| | 
* | commit a3f124c553931e32a02d2cdbdc2d690370d044b3
|/  Author: jamesengleback <jamesengleback@hotmail.co.uk>
|   Date:   Sat Jun 4 09:09:42 2022 +0100
|   
|       docs/ : nn.md - recomender systems start, index.md - gantt chart start, data-mining.md - small mod
|   
*   commit bc66710ccc01de1ef00f5a027a82ee641816e1de
|\  Merge: b5f5a4d 9dbf983
| | Author: jamesengleback <jamesengleback@hotmail.co.uk>
| | Date:   Thu Jun 2 17:55:49 2022 +0000
| | 
| |     Merge branch 'master' of https://github.com/jamesengleback/screening-fist
| |     merging onto gpu2
| |   
| *   commit 9dbf983e4f293e6d21758502a87cd6760e030f74
| |\  Merge: b4d49c7 32e1b8f
| | | Author: jamesengleback <james@engleback>
| | | Date:   Thu Jun 2 14:34:33 2022 +0000
...
```
