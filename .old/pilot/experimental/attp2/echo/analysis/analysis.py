import numpy as np
import pandas as pd
import plates

def main():
    df = pd.read_csv('compounds.csv') 
    smiles = df['SMILES'].to_list() + ['none', 'none']
    names = df['Item Name'].tolist() + ['none', 'none']


    concs = np.array([500 / 2**i for i in range(1,9)][::-1])
    
    blank = plates.UV384m3('data/blank.csv', 
                parser = 'ascii')
    wt = plates.UV384m3('data/3march-echo-wt-assay.csv', 
                parser = 'ascii')
    dm = plates.UV384m3('data/dm.csv', 
                parser = 'ascii')
    dm.report(save_dir = 'dm', 
                controlPlate = None, 
                concs = concs, 
                smiles = smiles, 
                names = names) 

    wt.report(save_dir = 'wt', 
                controlPlate = None, 
                concs = concs, 
                smiles = smiles, 
                names = names) 

if __name__ == '__main__':
    main()
