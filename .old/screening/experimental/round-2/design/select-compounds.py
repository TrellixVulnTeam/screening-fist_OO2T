import argparse
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters

def pick(smiles, n):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in mols]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    picker = SimDivFilters.MaxMinPicker()
    idx = picker.LazyPick(fn, len(smiles), n)
    return [smiles[i] for i in idx]

def lookup(smiles, df):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main(n):
    df = pd.read_csv(args.input)
    selection = lookup(pick(df['SMILES'], int(args.number)), df)
    selection.to_csv(args.output, 
            index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input')
    parser.add_argument('-o','--output')
    parser.add_argument('-n', '--number')
    args = parser.parse_args()
    main(args)
