#!/usr/bin/env python
import sys
import os
from tqdm import tqdm
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters

def get_mols(smiles):
    # duplicates
    usmiles = {i:Chem.MolFromSmiles(i) for i in set(smiles)}
    return [usmiles[i] for i in smiles]

def pick(smiles:list, n):
    picker = SimDivFilters.MaxMinPicker()
    #mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = {}
    for i in smiles:
        if i not in fps.keys():
            try:
                fps[i] = Chem.RDKFingerprint(Chem.MolFromSmiles(i))
            except:
                pass
    fps_ = list(fps.values())
    #fps = [Chem.RDKFingerprint(i) for i in mols.values()]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps_[i], fps_[j])
    if len(fps_) > n:
        idx = picker.LazyPick(fn, len(fps_), n)
        return [list(fps.keys())[i] for i in idx]

def main():
    path = '../data/pretraining/data/seq-smiles.csv'
    chunksize = 256
    n = 64
    df = pd.read_csv(path, 
                     chunksize=chunksize)
    smiles_sel = set()
    n_chunks = int(os.popen(f"wc -l {path}").read().split()[0]) // chunksize
    with tqdm(df, total=n_chunks) as bar:
        for chunk in bar:
            if bar is not None:
                sel = pick(chunk['smiles'].unique(), n)
                if sel is not None:
                    smiles_sel = smiles_sel.union(sel)
                    bar.set_postfix({'unique smiles':len(smiles_sel)})

    df = pd.read_csv(path)
    sys.stdout.write(f"{','.join(df.columns)}\n")
    for _ in range(8):
        for i in tqdm(smiles_sel):
            chunk = df.loc[df['smiles'] == i,:]
            if len(chunk) > 1:
                x = chunk.sample()
                if x is not None:
                    x.to_csv(sys.stdout, index=False, header=False)
            else:
                chunk.to_csv(sys.stdout, index=False, header=False)



if __name__ == '__main__':
    main()
