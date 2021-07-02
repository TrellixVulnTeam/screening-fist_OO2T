import pandas as pd
from picklist import SourcePlateCompound, Block, AssayPlate

def main():
    df = pd.read_csv('round1-cpds.csv')
    print(df)


if __name__ == '__main__':
    main()
