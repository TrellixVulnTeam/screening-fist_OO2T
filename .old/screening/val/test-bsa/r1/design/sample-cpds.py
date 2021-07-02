import pandas as pd
import cpds

def main():
    df = pd.read_csv('fda-lib.csv', index_col = 0)
    selected_smiles = cpds.cpds.pick(df.SMILES, 25)

    selected_df = pd.concat([df.loc[df.SMILES == i,:] for i in selected_smiles])

    selected_df.to_csv('selected_compounds.csv',index=False)

if __name__ == '__main__':
    main()
