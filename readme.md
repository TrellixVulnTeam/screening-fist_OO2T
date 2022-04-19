# Screening Fist

## Overview

**Screening Fist** is an enzyme compund screening operation where mutants of the cytochrome P450 BM3 are screened against many compounds for binding. 

## Aim

The aim of the operation is to provide a dataset of BM3 mutants, binding partners and their $K_d$ affinity values. The data is to be used by a pre-trained machine learning model to learn to predict binding affinity for a given BM3 mutant and a given compound. That model can then be used to virtually screen BM3 mutants for binding to a given compound by way of design.

## Protocol

The protocol is detailed in [protocols](docs/docs/protocol.md). Old work on assay development and pilot runs are in [.old](.old). 

#### Assay Overview

The assay detects binding between a purified P450 heme domain and a compound based on the shift in absorbance at 390 and 420 nm, analagous to traditional UV-Vis spectroscpy-based tutration experiments. No reactants are consumed and the assay gives a stable signal for > 4 hours with BM3 A82F/F87V. 

The assay is a plate assay where compounds are dispensed into the assay plate using a *Labcyte Echo* using picklists generated with [`echo`](https://github.com/jamesengleback/echo). 

Compounds can also be dispensed using a multichannel pipette from a serial dilution, but results are less accurate. 
P450 at around 10 uM in a suitable buffer is added to the plate and wells are read for absorbance between 250 and 800 nm. 
These wavelengths are used to closely replicate the traditional technique and also can be useful to debug issues.

Since some compounds absorb in the 250-800nm range, each compound concentration is paired with a matching concentration of that compound in buffer only. The absorbance of that compound in buffer is subtracted from the absorbance of that compound with protein, as with using a dual-beam spectrometer for P450 titrations.

In this work, the plate assay is 384 wells and is designed flexibly accordinng to the experimental requirements (e.g. het/miss screening with one high compound concentration or incramenting concentrations for determining a binding coefficient ($K_d$).

Data is processed with the package [`plates`](https://github.com/jamesengleback/plates) which contains tools for general (BMG) platereader data, and tools specific to this data.

For more information on the assay, look at [`protocols`](docs/docs/protocol.md).

## Screening Set

### BM3 mutants

BM3 mutants used in this work are selected from [1](wong) based on mutants with altered substrate scope. The advantage of this approach is that the screen is likely to have some hits and the disadvantage is that it is not an even distribution of sequences so care will have to be taken in training the model.

| ID | Mutations |
|----|-----------|
| wt |           |
| dm | A82F/F87V |

### Compounds 

The screening compound library is a [SelleckChem](https://www.selleckchem.com/) FDA approved screening library [L1300](https://www.selleckchem.com/screening/fda-approved-drug-library.html). Information and positions are in [lib](lib). 

> **note:** The screening library is out of date and has been used for other work in the past which has given me contamination and replicability concerns. This library was used because money was too tight to get a custom library, however work was done on library selection in [cpd-selection](dead)

## Repo Map

## References

[1](wong) 
