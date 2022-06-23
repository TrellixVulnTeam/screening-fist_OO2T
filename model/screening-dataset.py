#!/usr/bin/env python
import sys
import random
import pandas as pd

def smote(df, target='hit'):
    pos = df.loc[df[target],:]
    neg = df.loc[df[target] == False,:]
    if len(pos) > len(neg):
        x = neg.sample(len(pos), replace=True)
        o = pd.concat([neg, x])
    else:
        x = pos.sample(len(neg), replace=True)
        o = pd.concat([neg, x])
    return o.sample(frac=1).reset_index(drop=True)

def random_sample(df, frac=0.75):
    train = df.sample(frac=frac, replace=False)
    test = df.loc[[i for i in df.index if i not in train.index],:]
    return train, test

def main(args):
    for arg in args:
        if 'csv' in arg:
            df = pd.read_csv(arg)
        elif 'tsv' in arg:
            df = pd.read_csv(arg, delimiter='\t')
        else:
            raise Warning('bad file type maybe')
        data = smote(df)
        train, test = random_sample(data, 0.75)
        train.to_csv('data/screening-data.train.csv', index=False)
        test.to_csv('data/screening-data.test.csv', index=False)



if __name__ == '__main__':
    main(sys.argv[1:])
