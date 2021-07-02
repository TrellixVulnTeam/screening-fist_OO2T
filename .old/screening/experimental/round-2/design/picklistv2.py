from string import ascii_uppercase
import pandas as pd
import echo

def main():
    wells = [f'{i}{j}' for i in ascii_uppercase[:16] for  j in range(1,25)]
    df = pd.read_csv('round1-cpds.csv')

    # dmso
    dmso = echo.Cpd(name = 'dmso', vol = 10_000)
    dmso_src = echo.Src(ldv = False)
    for i in wells[:24]:
        dmso_src.fill(i, dmso, 60)
    
    # compounds
    cpds = [echo.Cpd(name = i, vol = 20) for i in df['Item Name']]
    cpd_src = echo.Src()
    for i, j in zip(wells, cpds):
        cpd_src.fill(i, j)

    # blocks and dest fill
    dest = echo.Dest()
    # block1 = echo.Block(cpds[0], dmso, dest.wells) # slice?


if __name__ == '__main__':
    main()
