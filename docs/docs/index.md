# Screening Fist

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
