#!/usr/bin/env python
import os
import argparse
from tqdm import tqdm

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from torch import Tensor
from einops import rearrange, repeat

from data import Data, DataTensors
from model import Model

def train(model, 
          data_loader, 
          epochs=8,
          lr=1e-6,
          #loss_fn='cross_entropy',
          proc=False,
          n_non_binders=0,
          **kwargs,
          ):
    #assert loss_fn in {'cross_entropy'}
    weight=Tensor([1.]*data_loader.batch_size \
                         + [0.]*data_loader.batch_size * n_non_binders)
    loss_fn = nn.BCELoss()#weight=weight)

    opt = torch.optim.Adam(model.parameters(), 
                           lr=lr,
                           )
    for epoch in range(epochs):
        # todo - data_loader sampler: importance sampling
        with tqdm(data_loader) as bar:
            for seq_, fingerprints_, hit_ in bar:
                seq = rearrange(seq_, 'b1 b2 l -> (b1 b2) l')
                fingerprints = rearrange(fingerprints_, 'b1 b2 l -> (b1 b2) l')
                hit = rearrange(hit_, 'b1 b2 -> (b1 b2)')
                yh = model(seq, fingerprints) # seq should be repeats, does it cache?
                if len(hit.shape) == 1:
                    y = rearrange(hit.float(), '(b c) -> b c', c=1)
                else:
                    y = hit.float()
                #print({'seq':seq.shape, 'fingerprints':fingerprints.shape, 'hit':hit.shape,
                #    'yh':yh.shape, 'y':y.shape})
                loss = loss_fn(yh, y)
                loss.backward()
                opt.step()
                opt.zero_grad()

                bar.set_postfix({'epoch':epoch,
                                 'loss':loss.detach().item(),
                                 })


def main(args):
    data = DataTensors(args.input, test=True, n_non_binders=3)
    #data = Data(args.input, test=True)
    a = len(data) // 4
    train_data, test_data = random_split(data, (len(data)-a, a))
    train_loader = DataLoader(train_data,
                              batch_size=2,
                              shuffle=True,
                              num_workers=0,
                              )
    test_loader = DataLoader(test_data,
                             batch_size=2,
                             shuffle=True,
                             num_workers=0,
                             )

    model = Model()

    train(model,
          train_loader,
          proc=False,
          n_non_binders=1,
          )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    args = parser.parse_args()
    assert args.input is not None
    assert os.path.isfile(args.input)
    main(args)
