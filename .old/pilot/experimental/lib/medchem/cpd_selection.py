import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit import SimDivFilters

def pick(smiles, n):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in smiles]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    picker = SimDivFilters.MaxMinPicker()
    idx = picker.LazyPick(fn, len(smiles), n)
    return pd.Series([smiles[i] for i in idx], index = idx)


def main():
    df = pd.read_csv('fda.csv') # smiles errors, scrape cas numbers
    print(pick(df['SMILES'], 10))

if __name__ == '__main__':
    main()
