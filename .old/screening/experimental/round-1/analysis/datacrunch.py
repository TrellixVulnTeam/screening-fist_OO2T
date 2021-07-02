import numpy as np
import pandas as pd
import plates

def main():
    concs = np.linspace(0,1,8)**3 * 500 # is k = 3?
    blank = plates.UV384m4('../data/blank-a.CSV')
    plate = plates.UV384m4('../data/a83f-a.CSV', 
            control = blank,
            concs = concs)
    for i in [1,2,3,4]:
        blk = plate.block(i, smiles = 'CCCCCC=O')
        #print(blk.df)
        #print(blk.norm)
        #print(blk.diff)
        #print(blk.response)
        print(blk.mm)
        #print(pd.DataFrame([plate.block(i).mm for i in range(1,24)]))
        #plates.report(blk)

if __name__ == '__main__':
    main()
