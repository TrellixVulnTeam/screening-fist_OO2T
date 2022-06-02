# Model

## Summary

$$
P_{binding} = fn(sequence, smiles)
$$

This page describes the deep learning model constructed for this project.
The model is designed to estimate the likelihood of a binding interaction between a given Cytochrome P450 sequence and a ligand SMILES code.
The intended end uses of the model are:

1. Virtually screening sequences for potential activity with a given compound.
2. Optimally design an enzyme-ligand screening experiment.

## Aim

- [ ] Build and train a model that: 
	- [x] Takes input enzyme sequence and a compound SMILES and outputs an estimate of binding likelihood interaction between the two.
- [ ] Pre-train prototype on large sequence-smiles dataset and save weights.
	- [x] Set up GPU machine
	- [ ] Run for $n$ iterations, save weights
	- [ ] Run again with validation
	- [ ] Outputs uncertainty


## Data

It was important to improve the sample efficiency of the model as far as possible, given the cost of lab data acquisition.
This was done by:

1. Application of a pre-trained model
2. Further pre-training the model on task-specific data

### Pre-Training Data

### Training Data

## Model 

The model aims to predict:
$$
P_{binding} = fn(sequence, smiles)
$$
Where $fn$ is a model that takes an input of a protein $sequence$ and a prospective ligands' SMILES code - $smiles$ and outputs $P_{binding}$ - an estimate of the probability that the given $sequence$ and $smiles$ bind to one another.

Given their prior success in chemical and protein sequence learning, a neural network model was chosen to build  $fn$.
The network can be split into three parts:


- **Sequence Embedding:** For a given protein $sequence$, output a tensor encoding a neural embedding $z_{sequence}$.
- **Chemical Embedding:** For a given chemical $smiles$ encoding, outputs an embedding $z_{smiles}$.
- **Prediction Head:** For the embeddings $z_{smiles}$ and $z_{sequence}$, output a prediction $P_{binding}$.

### Sequence Embedding
### Chemical Embedding
### Prediction Head

## Training

## Active Learning
