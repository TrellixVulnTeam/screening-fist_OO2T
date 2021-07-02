import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters


def filter_HL(smiles):
    # filter smiles by herbicide likeness rules Hao
    mw_cutoff_min = 100
    mw_cutoff_max = 435
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
    
    prop_filter = lambda s : props[s]['mw'] <= mw_cutoff_max \
                        and props[s]['mw'] >= mw_cutoff_min \
                        and props[s]['logp'] <= logp_cutoff \
                        and props[s]['hba'] <= hba_cutoff\
                        and props[s]['hbd'] <= hbd_cutoff\
                        and props[s]['rot'] <= rotatable_cutoff\
                        and props[s]['aroB'] <= n_aromatic_bonds_cutoff
    return [i for i in smiles if prop_filter(i)]

def lookup(smiles, df):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main():
    df = pd.read_csv('layouts.csv')
    herbicideLike = filter_HL(df['SMILES'])
    #herbicideLike = lookup(herbicideLike, df)
    print(len(herbicideLike))
    # herbicideLike.to_csv('herbicide-like.csv')

if __name__ == '__main__':
    main()
