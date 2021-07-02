import pandas as pd

def main():
    df = pd.read_excel('selleckchem-plate.xlsx',
            sheet_name = 'L1300-FDA-978cpds')
    df = df.loc[:,['Item Name','CatalogNumber', 'SMILES', 'Rack Number', 'Plate Location']]
    df.to_csv('layouts.csv')

if __name__ == '__main__':
    main()
