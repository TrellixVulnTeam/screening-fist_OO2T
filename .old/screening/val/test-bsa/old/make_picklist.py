import numpy as np
import pandas as pd
import echo

def make_block(wells, a, b, vol = 1.5, k = 3):
    # a nd b are compounds
    # vol is in µl
    # k is a skewing factor
    # return modified wells
    x = (np.arange(0,1, len(wells)) ** k) * vol
    for i,j in zip(wells, x):
         # change to xfer
         # required: parent plate/well, child plate/well
        i.fill(a, j)
        i.fill(b, vol - j)
    return wells


def main():
    df = pd.read_csv('selection.csv',index_col=0)
    dmso = echo.Cpd(name='dmso', vol = 1000)
    cpds = [echo.Cpd(name=i, vol=20) for i in df['CatalogNumber']]
    src_cpds = echo.Src()
    src_dmso = echo.Src(ldv=False)

    dest = echo.Dest()

    wells = make_block([dest.wells[i] for i in echo.hwells[:8]], cpds[0], dmso)
    print(wells)
    for i in wells:
        print(i.contents)

    # for i,j in zip(echo.hwells, cpds):
    #     src_cpds.fill(i, j, vol=12)

    # for i in echo.hwells[:5]:
    #     src_dmso.fill(i,dmso, 60)
    # 
    # for i in [1,2,3]:
    #     dest = echo.Dest()
    #     # make blocks


if __name__ == '__main__':
    main()
