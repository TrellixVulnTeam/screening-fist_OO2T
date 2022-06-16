#!/usr/bin/env python
import sys
import os
from hashlib import md5
from tqdm import tqdm
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from torch import Tensor, FloatTensor, cat, zeros
from model import Esm
from data import Data
from esm import Alphabet


class Data(Dataset):
    def __init__(self,
                 seq,
                 ):
        super().__init__()
        self.seq = seq
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, idx):
        return self.seq[idx]

TOKS=['<null_0>', '<pad>', '<eos>', '<unk>', 'L', 'A', 'G', 'V',
 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K',
 'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C',
 'X', 'B', 'U', 'Z', 'O', '.', '-', '<null_1>',
 '<cls>', '<mask>', '<sep>']

def main(args):
    hashfn = lambda s : md5(s.encode()).hexdigest()
    odir = 'embeddings'
    if not os.path.exists(odir):
        os.mkdir(odir)

    abc = Alphabet(TOKS)
    bc = abc.get_batch_converter()
    esm = Esm()
    esm.cuda()
    for arg in args:
        df = pd.read_csv(arg)
        # 800 is about all i can get away with on an RTX6000
        seq = [i for i in df['seq'].unique() if len(i) <= 800]
        data = Data(seq)
        loader = DataLoader(data, batch_size=16, num_workers=3)
        for batch in tqdm(loader):
            hashes = [hashfn(i) for i in batch]
            labels, sequences, tensors = bc(list(zip(hashes, batch)))
            tensors_ = cat([i.unsqueeze(0) for i in tensors])
            emb = esm(sequences).detach().cpu()
            for i,j in zip(emb, hashes):
                torch.save(i, os.path.join(odir, j +'.pt'))


if __name__ == '__main__':
    main(sys.argv[1:])
