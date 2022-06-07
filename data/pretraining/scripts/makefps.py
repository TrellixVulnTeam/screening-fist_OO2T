#!/usr/bin/env python
import sys
from multiprocessing.pool import ThreadPool

import pandas as pd
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.AllChem import RDKFingerprint

def smiles_fp(smiles):
    try:
        return list(RDKFingerprint(MolFromSmiles(smiles)))
    except:
        return []

def main(args):
    for arg in args:
        df = pd.read_csv(arg)
        smiles = df['smiles'].unique()
        with ThreadPool(128) as process_pool:
            fps_ = process_pool.map(smiles_fp, smiles)
        fps = pd.DataFrame(fps_)
        fps.index = smiles
        fps.to_csv(sys.stdout, header=False)



if __name__ == '__main__':
    main(sys.argv[1:])
