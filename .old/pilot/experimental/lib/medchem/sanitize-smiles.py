import time
from io import StringIO
import requests 
import pandas as pd
from tqdm import tqdm

def get_cid(cas):
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=pccompound&UID={cas}&rettype=SMILES'
    r = requests.get(url)
    print(r.status_code)
    if r.status_code == 200:
        print(r.text)
    if r.status_code == 414:
        print('url too long')

def search_name(name):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/CSV'
    r = requests.get(url)
    if r.status_code == 200:
        data = r.text
        df = pd.read_csv(StringIO(data))
        df['name'] = name
        return df
    else:
        print(name, r.status_code)
        pd.DataFrame([], columns = ['CID','CanonicalSMILES','name'])


def get_smiles(cid):
    pass

def main():
    df = pd.read_csv('fda.csv')
    cas = df['CAS Number']
    names = df['Product Name']


    results = pd.DataFrame([], columns = ['CID','CanonicalSMILES','name'])
    for i in tqdm(names):
        try:
            results = results.append(search_name(i))
        except :
            time.sleep(60)
            results = results.append(search_name(i))
    results.to_csv('fda-canonical-smiles.csv')

if __name__ == '__main__':
    main()
