#!/usr/bin/env python
import sys
import os
from hashlib import md5
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
        pad = lambda x, l : cat([x, zeros(l-len(x))])
        seqx = pad(Tensor(self.abc.encode(self.seq[idx])), 
                   self.max_seq_len).unsqueeze(0)
        #seqx = cat([seqx, zeros(self.max_seq_len - len(seqx))])
        fpx = smiles_fp(self.smiles[idx])
        if isinstance(self.hit[idx], (bool, str, int)):
            hitx = FloatTensor([self.hit[idx]])
        else:
            pass 
        seqfs = []
        fpfs = []
        for i in range(self.n_non_binders):
            i, j = random.randint(0, self.__len__()-1), \
                    random.randint(0, self.__len__()-1), 
            #seqf = pad(Tensor(self.abc.encode(self.seq[i])), 
            #           self.max_seq_len).unsqueeze(0)
            #seqfs.append(seqf)
            seqfs.append(seqx) # new - duplicating, lru_cache work?
            fpfs.append(smiles_fp(self.smiles[j]).unsqueeze(0))
        seqx = cat([seqx, *seqfs], dim=0)
        fpx = cat([fpx.unsqueeze(0), *fpfs], dim=0)
        hitx = cat([hitx, zeros(len(seqfs))], dim=0)
        #print({'seqx':seqx.shape, 'fpx':fpx.shape, 'hitx':hitx.shape})
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


class DataEmbeddings(Data):
    def __init__(self,
                 path,
                 embeddings_dir,
                 test=False,
                 cuda=True,
                 max_seq_len=512,
                 quiet=False,
                 **kwargs,
                 ):
        super().__init__(path=path, test=test, max_seq_len=max_seq_len, **kwargs)
        if cuda:
            self.embeds = {i.split('.')[0]:torch.load(os.path.join(embeddings_dir, i)) \
                            for i in os.listdir(embeddings_dir)}
        else:
            self.embeds = {i.split('.')[0]:torch.load(os.path.join(embeddings_dir, i),
                                                      map_location=torch.device('cpu')) \
                            for i in os.listdir(embeddings_dir)}
        self.max_seq_len = max_seq_len
        self.cuda = cuda
        self.hashfn = lambda s : md5(s.encode()).hexdigest()

        fails = []
        for i in set(self.seq):
            seq_hash = self.hashfn(i)
            if seq_hash not in self.embeds.keys():
                fails.append(i)
        if not quiet:
            if len(fails) > 0:
                raise Warning(f"Not all sequences embedded:\n{fails}")
        else:
            for i,j,k in zip(self.seq, self.smiles, self.hit):
                if i in fails:
                    self.seq.remove(i)
                    self.smiles.remove(j)
                    self.hit.remove(k)
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, idx):
        pad = lambda x, l : cat([x, zeros(self.max_seq_len+1-x.shape[0], x.shape[1])])
        seq_hash = self.hashfn(self.seq[idx])
        seqemb = self.embeds[seq_hash]
        seqx = pad(seqemb, self.max_seq_len)
        fpx = smiles_fp(self.smiles[idx])
        if isinstance(self.hit[idx], (bool, str, int)):
            hitx = FloatTensor([self.hit[idx]])
        else:
            pass 
        seqfs = []
        fpfs = []
        for i in range(self.n_non_binders):
            i, j = random.randint(0, self.__len__()-1), \
                    random.randint(0, self.__len__()-1), 
            #seqf = pad(Tensor(self.abc.encode(self.seq[i])), 
            #           self.max_seq_len).unsqueeze(0)
            #seqfs.append(seqf)
            seqfs.append(seqx) # new - duplicating, lru_cache work?
            fpfs.append(smiles_fp(self.smiles[j]).unsqueeze(0))
        seqx = cat([seqx, *seqfs], dim=0)
        fpx = cat([fpx.unsqueeze(0), *fpfs], dim=0)
        hitx = cat([hitx, zeros(len(seqfs))], dim=0)
        return seqx, fpx, hitx

@lru_cache(128)
def smiles_fp(smiles):
    return torch.Tensor(Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))).float()

def test(args):
    from torch.utils.data import DataLoader
    for arg in args:
        data = DataEmbeddings(arg, 
                              embeddings_dir='seq-smiles-embeds', 
                              test=False,
                              cuda=False,
                              quiet=True,
                              )
        print(len(data))
        data_loader = DataLoader(data,
                                 batch_size=32,
                                 )
        #data = DataTensors(arg, test=True)
        for seq, smiles, hit in tqdm(data_loader):
            #fp = smiles_fp(smiles)
            #print(seq, smiles, hit)
            #print(seq.shape, smiles.shape, hit.shape)
            pass

if __name__ == '__main__':
    test(sys.argv[1:])
