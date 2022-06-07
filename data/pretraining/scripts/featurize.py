#!/usr/bin/env python
import sys
from multiprocessing import pool                                                                   
from tqdm import tqdm
import pandas as pd
from rdkit.Chem import AllChem as Chem
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl

logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog('rdApp.error')

fns = ['CalcExactMolWt', 'CalcNumLipinskiHBD', 'CalcNumLipinskiHBA', 'CalcNumHBD', 'CalcNumHBA', 
 'CalcNumRotatableBonds', 'CalcNumRings', 'CalcNumAromaticRings', 'CalcNumSaturatedRings', 
 'CalcNumHeterocycles', 'CalcNumAromaticHeterocycles', 'CalcNumAromaticCarbocycles', 
 'CalcNumSaturatedCarbocycles', 'CalcNumAliphaticRings', 'CalcNumAliphaticHeterocycles', 
 'CalcNumHeteroatoms', 'CalcNumAmideBonds', 'CalcFractionCSP3']

def ft_smiles(smiles): 
    d = {} 
    try:
        m = Chem.MolFromSmiles(smiles)
        mh = Chem.AddHs(m)
        for i in fns:
            d[i] = Chem.__dict__[i](mh)
    except:
        pass
    return d 

 
def main():
    bs = 8
    fps = pd.read_csv('../data/seq-smiles.fp.csv', 
                      index_col=0, 
                      header=None, 
                      chunksize=bs,
                      )

    first = True
    for df in tqdm(fps):
        with pool.ThreadPool(bs) as ppool:
            ft = pd.DataFrame(list(ppool.map(ft_smiles, 
                                             df.index)),
                              index=df.index,
                              )
        if first:
            ft.to_csv('../data/seq-smiles-ft.csv')
            first = False
        else:
            ft.to_csv('../data/seq-smiles-ft.csv', mode='a', header=None)


if __name__ == '__main__':
    main()

