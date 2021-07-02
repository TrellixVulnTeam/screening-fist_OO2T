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
    selected_compounds = pd.read_csv('selection.csv', index_col=0).dropna()

    cpds = [echo.Cpd(i, vol=100) for i in selected_compounds.index]
    dmso = echo.Cpd('dmso', vol=1500)

    src1 = echo.SrcPlate(name='src1')
    # next free well: J1 (n=211)
    src_dmso = echo.SrcPlate(name='src_dmso', ldv=False)

    src_plate_fill = []
    layout = []
    
    # used up to inc f24
    #In [13]: echo.hwells.index('F24')
    #Out[13]: 143
    for i,j,k in zip(cpds, src1[144::2], src1[145::2]):
        j.fill(i.sample(12))
        k.fill(i.sample(12))
        src_plate_fill.append([j.loc,i.name])
        src_plate_fill.append([k.loc,i.name])
        layout.append({i.name:[j.loc,k.loc]})

    for i in src_dmso[211:211+19]:
        i.fill(dmso.sample(60))

    for plate in ['bsa-blank', 'bsa-test','ctrl-blank','ctrl-test']:
        dest = echo.DestPlate(name=f'dest-{plate}')
        for i, j in zip(range(0,24*8,8), cpds):
            xfer_block(dest[i:i+8], j,dmso)

    src_xfers = pd.DataFrame(src1.xfer_record)
    dmso_xfers = pd.DataFrame(src_dmso.xfer_record)
    src_plate_fill = pd.DataFrame(src_plate_fill)
    src_plate_fill.to_csv('src_plate_fill.csv',index=False)
    src_xfers.to_csv('cpd-src-picklist.csv')
    dmso_xfers.to_csv('dmso-src-picklist.csv')

    ### layout suff
    layout = pd.DataFrame(layout).melt().dropna()
    layout.columns=['CatalogNumber','Src Wells']
    layout.index=layout['CatalogNumber']
    fda = pd.read_csv('fda-lib-layout.csv', index_col=0)
    fda.index = fda['CatalogNumber']
    fda = fda.reindex(layout['CatalogNumber'])
    layout = pd.concat([fda, layout['Src Wells']], axis=1)
    layout = layout[['Item Name', 'CatalogNumber', 'Rack Number',  'Plate Location',  'Src Wells']]  
    layout.to_csv('src_layout.csv', index=False)



if __name__ == '__main__':
    main()
