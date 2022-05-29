#!/usr/bin/env python
import sys
import os

from rdkit import Chem
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from torch import cat, relu, sigmoid, Tensor, FloatTensor, LongTensor
from einops import rearrange

import esm
from data import Data

def fp(smiles):
    return torch.FloatTensor(\
            Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))\
            ).float().unsqueeze(0) # 0, 2048

class Esm(nn.Module):
    def __init__(self,
                 *,
                 model='esm1_t6_43M_UR50S',
                 ):
        super().__init__()
        pretrained = {'esm1_t6_43M_UR50S':esm.pretrained.esm1_t6_43M_UR50S}
        assert model in pretrained.keys()
        self.model, self.alphabet = pretrained[model]()
        self.batch_converter = self.alphabet.get_batch_converter()
    def __call__(self, seq):
        if isinstance(seq, str): # single prediction
            x = LongTensor(self.alphabet.encode(seq))
        elif isinstance(seq, (list, tuple)): # batch, flat list
            batch = [(i,j) for i,j in zip(range(len(seq)), seq)] # (id, seq), ...
            ids, seqs, x = self.batch_converter(batch)
        else:
            raise Warning(f"input types: str, list, tuple.\n{type(seq)}")
        z = self.model(x)
        return self.forward(x)
    def forward(self, x):
        return self.model(x)

class Fpnn(nn.Module):
    def __init__(self,
                 *,
                 fp_size=2048,
                 emb_size=32,
                 ):
        super().__init__()
        self.nn = nn.Sequential(nn.Linear(fp_size, emb_size),
                                nn.ReLU(),
                                )
    def __call__(self, smiles):
        if isinstance(smiles, str):
            fpt = fp(smiles)
        elif isinstance(smiles, Tensor):
            fpt = smiles
        elif isinstance(smiles, (tuple, list)):
            fpt = cat([fp(i) for i in smiles], dim=0)
        else:
            raise Warning(f"smiles input types: str, list, tuple.\n{type(smiles)}")
        print(fpt.shape)
        #return self.nn(fpt)

class Head(nn.Module):
    def __init__(self,
                 *,
                 emb_size=32,
                 ):
        super().__init__()
        self.nn = nn.Sequential(nn.Linear(emb_size, 1),
                                nn.ReLU(),
                                )
    def __call__(self, emb):
        yh = None
        return yh

class Model(nn.Module):
    def __init__(self,
                 *,
                 fpnn=Fpnn(),
                 head=Head(),
                 esm=Esm(),
                 ):
        super().__init__()
        self.fpnn = fpnn
        self.esm = esm.eval()
        self.head = head
    def __call__(self, seq, smiles):
        if isinstance(smiles, str):
            pass
        if isinstance(seq, str):
            pass
        fpz = self.fpnn(smiles)
        #seqz = self.esm(seq) # big, job killed
        #yh = self.head(seqz, fpz)
        #return yh

def main(arg='o.csv'):
    data = Data(arg, test=False)
    a = len(data) // 4
    train, test = random_split(data, (len(data)-a, a))
    train_loader = DataLoader(train,
                              batch_size=32,
                              shuffle=True,
                              num_workers=1,
                              )
    test_loader = DataLoader(test,
                             batch_size=32,
                             shuffle=True,
                             num_workers=1,
                             )

    esm = Esm()
    fpnn = Fpnn()
    head = Head()


if __name__ == '__main__':
    main()
