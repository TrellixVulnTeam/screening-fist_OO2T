# Data Mining

## Summary

It was necessary to pre-train the enzyme-ligand prediction model on a large and diverse dataset of enzyme sequences and their small molecule binding partners prior to training on the relatively small (~4000) dataset generated in this project.
Pre-training can equip a model with improved sample efficiency and generalizability, which are essential for using a sequence-ligand prediction model on expensive, lab-generated datasets.
A dataset of ~10<sup>6</sup> protein sequences and SMILES codes for their corresponding binding partners are gathered here from several sources.


## Aim

- Collect a large dataset of protein sequences and the SMILES codes for their binding partners.
- Collate into a dataset with sequence and chemical diversity.

## Background

### Sources

There are several sources of data containing protein sequences and information on their binding partners.

- [*Uniprot*](https://www.uniprot.org/) is a large, publicly available database containing protein sequences and a variety of functional annotations. 
*Uniprot* comprises three databases:
	- *UniProt Knowledgebase (UniprotKB)* contains protein sequences and their annotations.
	*UniprotKB* mainly comprises two collections:
		- *SProt:* - The SwissProt collection. Manually annotated and non-redundant protein sequences. 
		Annotations can include *Enzyme Comission (EC)* classification numbers, which can be linked to substrates and products. 
		At the time of collection (2022) the *SProt* collection contains 1,706,449 data points.
		- *Trembl:* - The TrEMBL collection. Automatic annotations, N data points.
	- *UniProt Archive (UniParc)* - an archive of sequences pooled from other public protein sequence databases.
	This collection does not contain annotations.
	- *UniProt Reference Clusters (UniRef)* - Three collections of clustered sequences from *UniprotKB* and *UniParc*. 
	Each collection is clustered with the CD-HIT algorithm and a different similarity cutoff: *UniRef100*, *UniRef90* and *UniRef50* use 100, 90 and 50% respectively.
	In the case of UniRef100, similarity measurements are made against sequence fragments.
	The clusters do not contain annotations.
Importantly: *Uniprot* files can be freely downloaded via an `ftp` server at [`https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/)

- [*KEGG*](https://www.genome.jp/kegg/) is another large source of proteins and biological function annotations, though download of the full dataset is restricted to paying customers.
*KEGG* does however have a *REST API* which can be used to query the database, though requests are handled somewhat slowly.
Despite this limitation, it was useful to cross-reference *EC* numbers against substrate and product identities.

- [*BindingDB*](https://www.bindingdb.org/bind/index.jsp) is a large collection of proteins and known small molecule binding partners mostly derived from high throughput drug screening datasets, which is reflected in the chemical diversity of the set.

- [*PDBBind*](http://pdbbind.org.cn/) contains protein sequences and structures and their binding ligands derived from the *Protein DataBank (PDB)*.

- [*BRENDA*](https://www.brenda-enzymes.org/) contains manual curations from primary literature for enzymes.

## Method

For this work, the *Uniprot SProt* and *BindingDB* collections were selected because they met the necessary requirements of containing protein sequences that can be matched to a chemical binding partner.
For data mining a cloud machine was hired (*Linode* `g6-dedicated-8`, 1TB Disk space, 8 `x86_64` cpus, 16 GB RAM, Debian 10, London).
Scripts used here are in `screening-fist/data/pretraining/scripts/`.

Large files and files in large numbers can prove slow to search, a problem that databases address.
For this work, a database was set up on the hired machine for storage and retrieval of the mined data.
*MongoDB*, an unstructured database, was chosen. 

*MongoDB* is unstructured in as far as it does not store tabular data, instead opting for a `json`-like schema which uses attribute:value arrays that themselves can contain more attribute:value arrays.
This offers considerable flexibility in data input from heterogeneous sources with varying presence, absence and categories of annotation.

*SProt* was gathered as a compressed `.dat` file using `screening-fist/data/pretraining/scripts/uniprot.sh`, which decompresses, parses entries and uploads them to the local *MongoDB* in a single *UNIX* pipe.
Using a single pipe saves the need to store the decompressed file on disk, which becomes significant in the case of the `TrEMBL` collection which has an uncompressed size of around 700GB as opposed to the 3.5GB of `sprot`.

*BRENDA* was used for its comprehensive list of *EC* numbers which would later be used for reference.
*BindingDB* was downloaded in full as a `.tsv` file.
Both the *BRENDA* and *BindingDB* sets were downloaded using `curl` commands imitating a web browser to bypass DDOS protection(?) and *BindingDB* was uploaded to the *MongoDB* instance.

The *EC* numbers from *BRENDA* were used with the *KEGG* REST API to collection compound identification numbers, which in turn were used to retrieve SMILES codes from PubChem using the script `ec.py`
`mkdataset1.py` was used to assemble a `csv` containing protein sequences and SMILES codes of their binding partners (1,587,674), which was further filtered to exclude SMILES codes and sequences that could not be parsed using RDKit and the tokenizer used in the sequence model.

!!! todo
	figure: hists of chemical properties, maybe a umap

Finally, an attempt to improve the chemical diversity of the pre-training dataset was made.
A diversity filter that aims to return a subset of a compound set by maximizing the pairwise Tanimoto similarity of the subset compound fingerprints was employed using the MaxMin algorithm.

!!! todo
	explain MaxMin

Unfortunately, MaxMin is $O(n^2)$ complexity, so is only feasible with relatively small batches.
64 compounds are selected this way from each batch of 512 from the filtered set, yielding `o3.csv` (1359834 data points and 908M uncompressed size).

`o3.csv` was compressed with `gzip` and loaded to an area of *Linode* object storage in Frankfurt, making it accessible to the *Linode* instances that would be used for model training.

`o3.csv` contained at least one invalid SMILES or sequences which had the costly effect of crashing the first model training run.
A filter for troublesome inputs was built and yeilded `o3f.csv` (x xMB).
and
For the purposes of model training and evaluation, the `o3f.csv` dataset was split 3:1 into a training and test datasets  `o3f.train.csv` and `o3f.test.csv` of sizes 907214 and 226804 respectively.
