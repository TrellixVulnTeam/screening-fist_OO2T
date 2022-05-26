#!/usr/bin/env python
import sys
import json
import pandas as pd

def proc(df):
    df = df.loc[df['Number of Protein Chains in Target (>1 implies a multichain complex)'] == 1,:]

    metrics = ['Ki (nM)',
               'IC50 (nM)',
               'Kd (nM)',
               'EC50 (nM)',
               #'kon (M-1-s-1)',
               #'koff (s-1)',
               ]
    smiles = df['Ligand SMILES']
    seq = df['BindingDB Target Chain  Sequence']

    o = []
    for i in metrics:
        x = df[i].astype(str).str.contains('[0-9\.]+').fillna(False) # hit
        gt = df[i].astype(str).str.startswith('>').fillna(False) #  miss
        o.append(x.astype(int) - gt.astype(int)) # 

    hits = pd.concat(o, axis=1).sum(axis=1).astype(bool)
    o = pd.concat([seq, smiles, hits], axis=1)
    o.columns = ['seq', 'smiles', 'hit']
    return o



def main(args):
    for arg in args:
        df = pd.read_csv(arg, delimiter='\t', on_bad_lines='skip', chunksize=1024)
        #df = pd.read_csv(arg, delimiter='\t', on_bad_lines='skip', nrows=10_000) # test
        for df_ in df:
            o = proc(df_).to_dict(orient='index')
            for i in o:
                json.dump(o[i], sys.stdout)


if __name__ == '__main__':
    main(sys.argv[1:])
