#!/usr/bin/env python
import sys
import random
from functools import lru_cache
from tqdm import tqdm
import pandas as pd
import torch
from torch.utils.data import Dataset
from torch import Tensor, FloatTensor, cat, zeros
from rdkit import Chem

from esm import Alphabet

TOKS=['<null_0>',
 '<pad>',
 '<eos>',
 '<unk>',
 'L',
 'A',
 'G',
 'V',
 'S',
 'E',
 'R',
 'T',
 'I',
 'D',
 'P',
 'K',
 'Q',
 'N',
 'F',
 'Y',
 'M',
 'H',
 'W',
 'C',
 'X',
 'B',
 'U',
 'Z',
 'O',
 '.',
 '-',
 '<null_1>',
 '<cls>',
 '<mask>',
 '<sep>']

class Data(Dataset):
    def __init__(self,
                 path,
                 test=False,
                 max_seq_len=512,
                 n_non_binders=0,
                 **kwargs,
                 ):
        super().__init__()
        self.path = path
        self.test = test
        self.n_non_binders = n_non_binders
        self.proc(max_seq_len=max_seq_len)
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, idx):
        if self.n_non_binders != 0:
            seqfs = []
            smilesfs = []
            for i in range(self.n_non_binders):
                i, j = random.randint(0, self.__len__()), \
                        random.randint(0, self.__len__()), 
                seqfs.append(self.seq[i])
                smilesfs.append(self.smiles[i])
            return (self.seq[idx], self.smiles[idx], self.hit[idx]),\
                   (seqfs, smilesfs, 0)
        else:
            return self.seq[idx], self.smiles[idx], self.hit[idx]
    def __repr__(self):
        return f"Dataset, {self.__len__()}"
    def proc(self, max_seq_len=None):
        if self.test:
            df = pd.read_csv(self.path, nrows=2048)
        else:
            df = pd.read_csv(self.path)
        if max_seq_len is not None:
            assert isinstance(max_seq_len, int)
            df = df.loc[df['seq'].str.len() <= max_seq_len, :]
        self.seq = list(df['seq'].str.upper())
        self.smiles = list(df['smiles'])
        self.hit = list(df['hit'].fillna(False).astype(int))
        assert len(self.seq) == len(self.smiles)
        assert len(self.seq) == len(self.hit)

class DataTensors(Data):
    def __init__(self,
                 path,
                 test=False,
                 max_seq_len=512,
                 **kwargs,
                 ):
        super().__init__(path=path, test=test, max_seq_len=max_seq_len, **kwargs)
        self.path = path
        self.test = test
        self.max_seq_len = max_seq_len
        self.proc(max_seq_len=max_seq_len)
        self.abc = Alphabet(TOKS)
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, idx):
        seqx = Tensor(self.abc.encode(self.seq[idx]))
        seqx = cat([seqx, zeros(self.max_seq_len - len(seqx))])
        fpx = smiles_fp(self.smiles[idx])
        if isinstance(self.hit[idx], (bool, str, int)):
            hitx = FloatTensor([self.hit[idx]])
        else:
            pass 
        if self.n_non_binders != 0:
            seqfs = []
            fpfs = []
            for i in range(self.n_non_binders):
                i, j = random.randint(0, self.__len__()), \
                        random.randint(0, self.__len__()), 
                seqf = Tensor(self.abc.encode(self.seq[i]))
                seqfs.append(cat([seqf, 
                                  zeros(self.max_seq_len - len(seqf))],
                                  dim=0))
                fpfs.append(smiles_fp(self.smiles[j]).unsqueeze(0))
            seqx = cat([seqx, *seqfs], dim=0)
            fpx = cat([fpx.unsqueeze(0), *fpfs], dim=0)
            hitx = cat([hitx, zeros(len(seqfs))], dim=0)
        return seqx, fpx, hitx
    def __repr__(self):
        return f"Dataset, {self.__len__()}"
    #def proc(self, max_seq_len=None):
    #    if self.test:
    #        df = pd.read_csv(self.path, nrows=2048)
    #    else:
    #        df = pd.read_csv(self.path)
    #    if max_seq_len is not None:
    #        assert isinstance(max_seq_len, int)
    #        df = df.loc[df['seq'].str.len() <= max_seq_len, :]
    #    self.seq = list(df['seq'])
    #    self.smiles = list(df['smiles'])
    #    self.hit = list(df['hit'].fillna(False).astype(int))
    #    assert len(self.seq) == len(self.smiles)


@lru_cache(128)
def smiles_fp(smiles):
    return torch.Tensor(Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))).float()

def main(args):
    from torch.utils.data import DataLoader
    for arg in args:
        data = DataTensors(arg, test=False)
        data_loader = DataLoader(data,
                                 batch_size=32,
                                 )
        #data = DataTensors(arg, test=True)
        for seq, smiles, hit in tqdm(data_loader):
            #fp = smiles_fp(smiles)
            #print(seq, smiles, hit)
            pass

if __name__ == '__main__':
    main(sys.argv[1:])
