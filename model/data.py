#!/usr/bin/env python
import sys
from functools import lru_cache
import pandas as pd
import torch
from torch.utils.data import Dataset
from rdkit import Chem


class Data(Dataset):
    def __init__(self,
                 path,
                 test=False,
                 ):
        super().__init__()
        self.path = path
        self.test = test
        self.proc()
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, idx):
        return self.seq[idx], self.smiles[idx], self.hit[idx]
    def __repr__(self):
        return f"Dataset, {self.__len__()}"
    def proc(self):
        if self.test:
            df = pd.read_csv(self.path, nrows=2048)
        else:
            df = pd.read_csv(self.path)
        self.seq = list(df['seq'])
        self.smiles = list(df['smiles'])
        self.hit = list(df['hit'])
        assert len(self.seq) == len(self.smiles)

@lru_cache(128)
def smiles_fp(smiles):
    return torch.Tensor(Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))).float()

def main(args):
    for arg in args:
        data = Data(arg, test=False)
        for seq, smiles in data:
            fp = smiles_fp(smiles)

if __name__ == '__main__':
    main(sys.argv[1:])
