import heapq
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters
from tqdm import tqdm

def filter_HL(smiles):
    # filter smiles by herbicide likeness rules Hao
    # todo: n aromatic bonds < 17
    mw_cutoff = 435
    logp_cutoff = 6
    hba_cutoff = 6
    hbd_cutoff = 2
    rotatable_cutoff =  9
    n_aromatic_bonds_cutoff = 17
    
    makemols = lambda smiles : Chem.AddHs(Chem.MolFromSmiles(smiles))
    n_aromatic_bonds = lambda m : sum([b.GetIsAromatic() for b in m.GetBonds()])
    mols = [makemols(i) for i in smiles]
    props = {s:{'mw':Chem.CalcExactMolWt(i), 
                'logp': Crippen.MolLogP(i),
                'hba': Chem.CalcNumLipinskiHBA(i),
                'hbd': Chem.CalcNumLipinskiHBD(i),
                'rot': Chem.CalcNumRotatableBonds(i),
                'aroB': n_aromatic_bonds(i)}
                for i,s in zip(mols, smiles)}
    
    prop_filter = lambda s : props[s]['mw'] <= mw_cutoff \
                        and props[s]['logp'] <= logp_cutoff \
                        and props[s]['hba'] <= hba_cutoff\
                        and props[s]['hbd'] <= hbd_cutoff\
                        and props[s]['rot'] <= rotatable_cutoff\
                        and props[s]['aroB'] <= n_aromatic_bonds_cutoff

    return [i for i in smiles if prop_filter(i)]


def pick(smiles, n):
    picker = SimDivFilters.MaxMinPicker()
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in mols]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    idx = picker.LazyPick(fn, len(smiles), n)
    return [smiles[i] for i in idx]


def pick_set(filteredSMILES, plates, n):
    selectionSMILES = pick(filteredSMILES, n)
    selection = pd.concat([plates.loc[plates['CanonicalSMILES'] == i, :] for i in selectionSMILES])
    locations = selection[['Plate','Well','Product Name']]
    return locations.sort_values(['Plate','Well'])

def main():
    smiles = pd.read_csv('lib/fda-canonicalSMILES-deduplicated.csv', index_col=-1)['CanonicalSMILES']
    

    plates = pd.read_csv('lib/fda.csv')
    in_stock = ['HY-L022-1', 'HY-L022-2']
    lib = plates.loc[plates['Plate'].isin(in_stock), :]

    canonicalSMILES = []
    for i in lib['Product Name']:
        if i not in smiles.index:
            canonicalSMILES.append(None)
        else:
            canonicalSMILES.append(smiles.loc[i])

    lib = pd.concat([lib, pd.Series(canonicalSMILES, name = 'CanonicalSMILES')], axis = 1, join = 'inner')
    lib.drop(([i for i, j in zip(lib.index, lib['CanonicalSMILES']) if j is None]), inplace = True)


    filtered = filter_HL(lib['CanonicalSMILES'])
    selection = pick_set(filtered, lib, 48)
    print(selection)
    selection.to_csv('selection.csv', index=False)

if __name__ == '__main__':
    main()
