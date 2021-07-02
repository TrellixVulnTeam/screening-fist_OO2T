import pandas as pd
import cpds

def get_by_smiles(df, smiles):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main():
    df = pd.read_csv('fda-lib-layout.csv', index_col=0)
    smiles = df['SMILES']
    selection_smiles = cpds.pick(smiles, 24)

    selection = get_by_smiles(df, selection_smiles)
    selection.to_csv('selection.csv')

if __name__ == '__main__':
    main()
