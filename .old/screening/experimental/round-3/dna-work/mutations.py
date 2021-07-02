import re
from tqdm import tqdm 
import pandas as pd 

import mxn
from bm3 import bm3, orf


COMPLIMENT = {'A':'T','T':'A','C':'G','G':'C'}

reverse_comp = lambda s : ''.join([COMPLIMENT[i.upper()] for i in s[::-1]])

def parse_mutation(m):
    return int(re.findall('\d+', m)[0]), re.findall('\w', m)[-1] 

def main():
    with open('mutations.txt', 'r') as f:
        mutations = f.read().split('\n')[:-1]

    output = []
    for i in tqdm(mutations):
        cds = mxn.CDS(orf, bm3)
        try:
            pos, aa = parse_mutation(i)
            cds.mutate(pos, aa)
            output.append(pd.DataFrame(cds.primers).T)
        except:
            print(i)
    df = pd.concat(output)
    df.to_csv('primers.csv')





if __name__ == '__main__':
    main()
