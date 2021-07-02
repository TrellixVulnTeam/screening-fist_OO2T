import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

def main():
    df = pd.read_csv('fda-canonical-smiles.csv')[['name','CanonicalSMILES']].drop_duplicates()

    repeats = df['name'].value_counts().loc[df['name'].value_counts() > 1]

    mols = []
    names = []
    for i in repeats.index:
        x = df.loc[df['name'] == i, :]
        for j in x['CanonicalSMILES']:
            mols.append(Chem.MolFromSmiles(j))
            names.append(i)
    im = Draw.MolsToGridImage(mols, legends = names)
    im.save('duplicated-smiles.png')

if __name__ == '__main__':
    main()
