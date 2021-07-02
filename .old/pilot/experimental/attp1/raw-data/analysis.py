import os
import numpy as np
import pandas as pd
import plates

def main():
    # smiles and names
    df = pd.read_csv('compounds.csv', index_col=0)
    d12_names = df['Item Name'][:24]
    d12_smiles = df['SMILES'][:24]
    d34_names = df['Item Name'][25:49]
    d34_smiles = df['SMILES'][25:49]
    concs = np.array([500 / 2**i for i in range(1,9)][::-1])

    # reports
    dfs = []
    for i, j in zip(['wt-d1.CSV','a82-d1,2.CSV','dm-d1,2.CSV'], ['wt','a82f','a82f-f87v']):
        plate = plates.UV384m1(i,
                name = j + 'd12', 
                concs = concs, 
                smiles = d12_smiles, 
                names = d12_names)
        dfs.append(plates.reportPlate(plate, save_dir = j + 'd12'))

    for i, j in zip(['wt-d3+4.CSV','a82-d3,4.CSV','dm-d3,4.CSV'], ['wt','a82f','a82f-f87v']):
        plate = plates.UV384m1(i,
                name = j + '34', 
                concs = concs, 
                smiles = d34_smiles, 
                names = d34_names)
        dfs.append(plates.reportPlate(plate, save_dir = j + 'd34'))

    df = pd.concat(dfs).reset_index(drop=True)
    df.to_csv('screeningSummary.csv')

    # anomaly detection
    # csv aggregation

if __name__ == '__main__':
    main()
