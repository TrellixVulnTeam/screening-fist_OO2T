import numpy as np
import pandas as pd
import echo


def square_layout():
    # 
    rows = [echo.hwells[i:i+24] for i in range(0,384,24)]
    x=[]
    for i,j in zip(rows[::2],rows[1::2]):
        for k in range(0,24,2):
            x.append(i[k:k+2] + j[k:k+2])
    return x

def main():
    # plate 1 a1:a1
    # plate 2 a1:b1
    # plate 3 a1:a2
    # plate 4 a1:b2

    compounds = pd.read_csv('fda.csv', index_col=0)

    cpds = [echo.Cpd(i, vol=100) for i in compounds['CatalogNumber']]
    cpds_blocks_of_384 = [cpds[i:i+384] for i in range(0, len(cpds), 384)]
    dmso = echo.Cpd('dmso', vol=1500)

    src_dmso = echo.SrcPlate(name='src_dmso', ldv=False)


    src_plates = []
    for i,j in enumerate(cpds_blocks_of_384,1):
        src = echo.SrcPlate(name=f'src{i}')
        for k,l in zip(j, square_layout()):
            for m in l:
                src[m].fill(k.sample(12))
        src_plates.append(src)
    

    for i in src_dmso[192:192+19]:
        i.fill(dmso.sample(60))

    for plate in [1,2,3,4]:
        dest = echo.DestPlate(name=f'dest-{plate}')
        for i,j in zip(dest, cpds):
            j.xfer(12,i)


    #src_xfers = pd.DataFrame(src1.xfer_record)
    #dmso_xfers = pd.DataFrame(src_dmso.xfer_record)
    #src_plate_fill = pd.DataFrame(src_plate_fill)
    #src_plate_fill.to_csv('src_plate_fill.csv',index=False)
    #src_xfers.to_csv('cpd-src-picklist.csv')
    #dmso_xfers.to_csv('dmso-src-picklist.csv')


if __name__ == '__main__':
    main()
