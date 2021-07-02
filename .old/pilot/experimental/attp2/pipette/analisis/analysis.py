import numpy as np
import pandas as pd
import plates

def main():
    layout = pd.read_csv('compounds.csv', index_col=0)
    smiles = layout['SMILES']
    names = layout['Item Name']
    concs = np.array([500 / 2**i for i in range(1,9)])
    blank_data = plates.UV384m2('data/blank-pipette.CSV')
    wt_data = plates.UV384m2('data/wt-pipette.CSV') 
    dm_data = plates.UV384m2('data/dm-pipette.csv') 
    wt_data.report(smiles = smiles, names = names, concs = concs)
    dm_data.report(smiles = smiles, names = names, concs = concs)

if __name__ == '__main__':
    main()
