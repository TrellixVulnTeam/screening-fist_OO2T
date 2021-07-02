import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def process(df):
    df.index = df.iloc[:,1]
    df = df.iloc[:,2:]
    return df

def plot(df, name = None):
    plt.figure(figsize=(15,5))
    sns.heatmap(np.random.randn(8,12), 
            annot=df, 
            fmt = '',
            cbar = False,
            annot_kws={'fontsize':10,
                'rotation':30})
    plt.xticks(range(1,13), range(1,13))
    plt.yticks(range(8), list('abcdefgh'), rotation = 0)
    plt.title(name)
    if name != None:
        plt.savefig(name + '.png')
    else:
        plt.show()



def main():
    csvs = [i for i in os.listdir() if 'hy-' in i and 'csv' in i]
    dfs = [pd.read_csv(i, header=None) for i in csvs]
    dfs = [process(i) for i in dfs]
    for i, j in tqdm(zip(dfs, csvs)):
        try:
            plot(i, name = j.split('.')[0])
        except  Exception as e:
            print(e)

if __name__ == '__main__':
    main()
