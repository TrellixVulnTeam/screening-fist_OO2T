# Enzymes

## Summary

This page contains the methods for producing the enzymes used in this screening program.
The enzymes are variants of the Cytochrome P450 BM3:

| ID | Mutations | PDB  |
|----|-----------|------|
| WT |           | 1BU7 |
|A82F|   A82F    | 4KEW |
| DM | A82F/F87V | 4KEY |
|1YQO|   T268A   | 1YQO |
|1YQP|   T268N   | 1YQP |

The page shows the method used to create the mutant BM3 expression plasmid DNA, 
expression of the mutants in *E. coli* and their purification.

## Aims

- Create expression plasmids containing the target mutants from an in-house starting point - [bm3-wt.gb]().
	- Sequence the plasmids to confirm they carry the mutations
- Express the mutants in *E. coli* using those plasmids.
- Purify the mutant protein from the *E. coli* harvest.


## DNA

### Starting Material

An heirloom BM3 Wild-type (heme domain) expression plasmid, [bm3-wt.gb](), 
was inherited and used as the basis for DNA work in this project.
The plasmid is a **?** pET15(?) expression vector where the BM3 gene has a 6xHis purification tag at the N-terminus,
flanked by a T7 promoter and terminator which leads to high yields in strains of *E. coli* containing the T7 RNA polymerase.
The plasmid also encodes ampicillin resistance and a ? replication origin which leads to a low copy number.

!!! todo
	get plasmid info

### Primer design and Acquisition

Mutations were introduced to the wild-type sequence via Polymerase Chain Reaction (PCR)-based site-directed mutagenesis.
Two methods were considered for this task based on commercially available kits, where each imposes different constraints on primer design.
Efforts were made to automate primer design as far as possible with scalability in mind.

The PCR kits used were:

1. *New England Biolabs (NEB) Q5 mutagenesis kit* - which requires that primers facilitate cloning of a linear DNA strand from the circular template plasmid and mutation payloads are carried in the tail of one primer.
The kit includes a cocktail of the DNAse *DPN1*, which disassembles template plasmid methylated in *E. coli* and a kinase and ligase that work to join the ends of the linear DNA into a circular plasmid.
The reaction is restricted to one payload.

2. *Agilent Quickchange mutagenesis kit* - which requires a pair of overlapping primers that carry the mutation payload in the mid-section.
This cloning method produces circular DNA carrying the targeted changes. 
It has the advantage of allowing multiple payloads carried by multiple primer sets.

Two important considerations based on the template sequence are:

1. Adenine-Thymine (AT) richness of the template sequence. Compared to cytosine and guanine (C and G), A and T bind to their complimentary bases weakly.
This results in weak primer binding to the template sequence, measurable by a low primer *melting temperature* $T_m$.
To compensate, primers must be longer than they otherwise would be for a sequence richer in CG, which increases their cost and their chance of self-binding.
The template sequence used here is AT-rich - at $x%$
2. Repetitions and palindromic regions of the template sequence. 
If the sequence surrounding a mutation target area contains these features, then the likelihood of *mis-priming* by binding to an off-target sequence area is high, so too is the likelihood of a non-functional, self-binding primer.

!!! todo 
	- get AT and repetition content of template sequence
	- talk about auto primer design method


### PCR and Work Up

#### Materials

- Primers were ordered from *Eurofins Genomics* and diluted to a stock concentration of 10 uM??
- pcr kit
- starting plasmid
- dh5a 

#### Method



### Sequencing

Purified plasmid DNA ostensibly conataining the target mutations, having been harvested and purified from DH5a *E. coli* cells, was shipped to *Eurofins Genomics* for sequencing using their *TubeSeq* service, which uses a variant of Sanger Sequencing.
Sequencing primers for this matched the T7 promoter and terminator and provided coverage of the targetted region.

!!! todo
	sequencing analysis

## Expression

Having been sequenced and confirmed to carry the target mutations, the mutant plasmids were used to produce the mutant protein *en masse* via a *BL21 DE3 E. coli* strain, which contains a T7 RNA polymerase under the control of a *lac* promoter.

#### Materials

- Expression plasmid encoding mutant P450 BM3
- *BL21 DE3 E. coli* - NEB. This domesticated *E.coli* strain is shipped in a transformation buffer.
- Auto-induction *Terrific Broth* (TB) media, which contains glucose and a lactose analog. 
The lactose analog triggers expression of T7 RNA-polymerase in *BL21 DE3 E. coli* and the subsequent expression of the target protein between the T7 promoter and terminator regions. 
The glucose inhibits this until it is consumed by the cells, which allows them to multiply to sufficient numbers before diverting energy to production of the target protein.
- Ampicillin - the antibiotic for which resistance is encoded in the target plasmid, ensuring that all cells in the growth media contain this resistance.
Assuming no ampicillin-resistant contaminants, all cells should be *BL21 DE3 E. coli* containing the target plasmid.
- $\Delta$ amino-levulnic acid (D-ALA) - a precursor to heme, ensuring heme availability for the large amount of BM3.

#### Method

## Purification
