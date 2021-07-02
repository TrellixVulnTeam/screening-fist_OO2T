import numpy as np
import pandas as pd
import echo

# cpd needs to know which wells are its children 
# and send xfer command
# ? where to register xfer ?

def xfer_block(wells, a, b, vol = 1.5, k = 3):
    # wells: iterable of echo.Wells
    # a, b echo.Cpd
    # vol µl
    # k skew factor
    vols_a = (np.linspace(0,1,len(wells)) ** k) * vol
    vols_b = vol - vols_a
    for i,j,k in zip(wells, vols_a, vols_b):
        a.xfer(i,j)
        b.xfer(i,k)



def main():
    selected_compounds = pd.read_csv('selected_compounds.csv', index_col=0)

    cpds = [echo.Cpd(i, vol=100) for i in selected_compounds['CatalogNumber']]
    dmso = echo.Cpd('dmso', vol=1500)

    src1 = echo.SrcPlate(name='src1')
    # next free well: I1 (n=192)
    src_dmso = echo.SrcPlate(name='src_dmso', ldv=False)

    src_plate_fill = []
    
    for i,j,k in zip(cpds, src1[::2], src1[1::2]):
        j.fill(i.sample(12))
        k.fill(i.sample(12))
        src_plate_fill.append([j.loc,i.name])
        src_plate_fill.append([k.loc,i.name])

    for i in src_dmso[192:192+19]:
        i.fill(dmso.sample(60))

    for plate in ['bsa-blank', 'bsa-test','ctrl-blank','ctrl-test']:
        dest = echo.DestPlate(name=f'dest-{plate}')
        for i, j in zip(range(0,24*8,8), cpds):
            xfer_block(dest[i:i+8], j,dmso)
        dest.map.to_csv(f'{plate}-map.csv')

    src_xfers = src1.xfer_record
    dmso_xfers = src_dmso.xfer_record
    src_xfers.to_csv('cpd-src-picklist.csv')
    dmso_xfers.to_csv('dmso-src-picklist.csv')


if __name__ == '__main__':
    main()
