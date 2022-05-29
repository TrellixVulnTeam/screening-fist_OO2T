#!/usr/bin/env python
import os
import argparse
from tqdm import tqdm

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from data import Data
from model import Model

def train(model, 
          data_loader, 
          epochs=8,
          lr=1e-6,
          loss_fn='cross_entropy',
          **kwargs,
          ):
    assert loss_fn in {'cross_entropy'}
    if loss_fn == 'cross_entropy':
        loss_fn = nn.CrossEntropyLoss()

    opt = torch.optim.Adam(model.parameters(), 
                           lr=lr,
                           )
    for epoch in range(epochs):
        # todo - data_loader sampler: importance sampling
        with tqdm(data_loader) as bar:
            for seq, smiles, hit in bar:
                yh = model(seq, smiles)
                loss = loss_fn(yh, hit)
                loss.backward()
                opt.step()
                opt.zero_grad()

                bar.set_postfix({'epoch':epoch,
                                 'loss':loss.detach().item(),
                                 })


def main(args):
    data = Data(args.input, test=True)
    a = len(data) // 4
    train_data, test_data = random_split(data, (len(data)-a, a))
    train_loader = DataLoader(train_data,
                              batch_size=32,
                              shuffle=True,
                              num_workers=1,
                              )
    test_loader = DataLoader(test_data,
                             batch_size=32,
                             shuffle=True,
                             num_workers=1,
                             )

    model = Model()

    train(model,
          train_loader,
          )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    args = parser.parse_args()
    assert args.input is not None
    assert os.path.isfile(args.input)
    main(args)
