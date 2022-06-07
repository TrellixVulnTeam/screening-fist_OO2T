#!/usr/bin/env python
import sys
import os
from functools import lru_cache

from rdkit import Chem
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from torch import cat, relu, sigmoid, Tensor, FloatTensor, LongTensor
from einops import rearrange
from einops.layers.torch import Rearrange

import esm
from data import Data, DataTensors

def fp(smiles):
    return torch.FloatTensor(\
            Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))\
            ).float().unsqueeze(0) # 0, 2048

class Skip(nn.Module):
    def __init__(self, 
                 size,
                 *args):
        super().__init__()
        self.nn = nn.Sequential(nn.Linear(size,size),
                                nn.Dropout(0.2),
                                nn.ReLU(),
                                nn.BatchNorm1d(size),
                                )
        self.bn = nn.BatchNorm1d(size)
    def forward(self, x):
        return self.bn(self.nn(x) + x)

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
        device = next(self.model.parameters()).device
        if isinstance(seq, str): # single prediction
            x = LongTensor(self.alphabet.encode(seq))
        elif isinstance(seq, (list, tuple)): # batch, flat list
            batch = [(i,j) for i,j in zip(range(len(seq)), seq)] # (id, seq), ...
            ids, seqs, x = self.batch_converter(batch)
        elif isinstance(seq, (Tensor, FloatTensor, LongTensor)):
            x = seq.int()
        else:
            raise Warning(f"input types: str, list, tuple.\n{type(seq)}")
        x = x.to(device)
        return self.forward(x)
    def forward(self, x):
        return self.model(x)['logits'] # size : b l d=35

class SeqPool(nn.Module):
    def __init__(self,
                 input_channels=35,
                 conv_channels=35,
                 num_conv_layers=3,
                 kernel_size=9,
                 stride=3,
                 num_lstm_layers=2,
                 lstm_hs=32,
                 ):
        super().__init__()
        self.conv_channels = conv_channels
        self.num_conv_layers = num_conv_layers
        self.kernel_size = kernel_size
        self.stride = stride
        self.num_lstm_layers = num_lstm_layers
        self.lstm_hs = lstm_hs

        self.nn = nn.Sequential(nn.Conv1d(in_channels=input_channels,
                                          out_channels=conv_channels,
                                          kernel_size=1,
                                          stride=1),
                                nn.ReLU(),
                *[nn.Sequential(nn.Conv1d(in_channels=conv_channels,
                                          out_channels=conv_channels,
                                          kernel_size=kernel_size,
                                          stride=stride),
                                nn.ReLU(),)
                                for i in range(num_conv_layers)],
                Rearrange('b d l -> b l d'),
                nn.LSTM(input_size=conv_channels,
                    hidden_size=lstm_hs,
                    num_layers=num_lstm_layers,
                    batch_first=True,
                    ),
                )
    def __call__(self, seqz):
        zh = self.forward(seqz)
        return zh
    def forward(self, z):
        z = rearrange(z, 'b l d -> b d l')
        output, (hn, cn) = self.nn(z)

        # hn shape: 
        zh = rearrange(hn, 'l b d -> b (l d)')
        #print({'hn':hn.shape, 'zh':zh.shape})
        return zh

class Fpnn(nn.Module):
    def __init__(self,
                 *,
                 fp_size=2048,
                 emb_size=32,
                 n_layers=3,
                 ):
        super().__init__()
        self.nn = nn.Sequential(nn.Linear(fp_size, emb_size),
                                nn.Dropout(0.2),
                                nn.ReLU(),
                                nn.BatchNorm1d(emb_size),
                                *[Skip(emb_size) for _ in range(n_layers)],
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
        return self.nn(fpt)

class Head(nn.Module):
    def __init__(self,
                 *,
                 emb_size=96,
                 n_layers=3,
                 ):
        super().__init__()
        self.nn = nn.Sequential(\
                *[Skip(emb_size) for _ in range(n_layers)],
                #*[nn.Sequential(nn.Linear(emb_size, emb_size),
                #                nn.ReLU()) 
                #        for _ in range(n_layers)],
                nn.Linear(emb_size, 1),
                nn.Sigmoid(),
                                )
    def __call__(self, seqz, fpz):
        #print({'seqz':seqz.shape, 'fpz':fpz.shape})
        #seqzz = rearrange(seqz, 'b l d -> b (l d)')
        z = cat([seqz, fpz], dim=1)
        return self.forward(z)
    def forward(self, z):
        return self.nn(z)

class Model(nn.Module):
    def __init__(self,
                 *,
                 esm=Esm(model='esm1_t6_43M_UR50S'),
                 seqpool=SeqPool(input_channels=35,
                                 conv_channels=35,  
                                 num_conv_layers=3,
                                 kernel_size=9,
                                 stride=3,
                                 num_lstm_layers=2,
                                 lstm_hs=32,
                                 ),
                 fpnn=Fpnn(fp_size=2048,
                           emb_size=32,
                           ),
                 head=Head(emb_size=96,
                           n_layers=3,
                           ),
                 ):
        super().__init__()
        self.esm = esm.eval()
        self.seqpool = seqpool
        self.fpnn = fpnn
        self.head = head
    def __call__(self, seq, smiles):
        if isinstance(smiles, str):
            pass
        if isinstance(seq, str):
            pass
        fpz = self.fpnn(smiles)
        seqz = self.esm(seq) 
        seqzz = self.seqpool(seqz)
        #print({'seqz':seqz.shape,'seqzz':seqzz.shape, 'fpz':fpz.shape})
        yh = self.head(seqzz, fpz)
        return yh

def main(arg='o.csv'):
    data = DataTensors(arg, test=False)
    a = len(data) // 4
    train, test = random_split(data, (len(data)-a, a))
    train_loader = DataLoader(train,
                              batch_size=8,
                              shuffle=True,
                              num_workers=1,
                              )
    test_loader = DataLoader(test,
                             batch_size=8,
                             shuffle=True,
                             num_workers=1,
                             )

    esm = Esm()
    fpnn = Fpnn()
    head = Head()


if __name__ == '__main__':
    main()
