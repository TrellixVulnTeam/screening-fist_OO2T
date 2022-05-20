#!/usr/bin/env python
import os
import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import Draw

def main():
    df = pd.read_csv('layouts.csv', index_col=0)
    names = df['Item Name']
    ids = df['CatalogNumber']
    smiles = df['SMILES']

    if not os.path.exists('img'):
        os.mkdir('img')

    for i,j,k in tqdm(zip(names, ids, smiles), total=len(df)):
        mol = Chem.MolFromSmiles(k)
        im = Draw.MolToImage(mol)
        im.save(os.path.join('img', j+'.png'))



if __name__ == '__main__':
    main()
