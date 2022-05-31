#!/usr/bin/env python
import os
import argparse
from tqdm import tqdm

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from torch import Tensor
from einops import rearrange, repeat
import wandb

from data import Data, DataTensors
from model import Model, Esm, Head, Fpnn, SeqPool

def train(model, 
          data_loader, 
          epochs=8,
          lr=1e-6,
          #loss_fn='cross_entropy',
          proc=False,
          n_non_binders=0,
          save_path=None,
          cuda=False,
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
                if cuda:
                    seq_, fingerprints_, hit_ = seq_.cuda(), fingerprints_.cuda(), hit_.cuda()
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
                wandb.log({'epoch':epoch,
                          'loss':loss.detach().cpu().item(),
                           })
            if save_path is not None:
                torch.save(model.state_dict(), 
                           os.path.join(save_path, 
   f"{save_path.split('/')[-1]}_e{epoch}l{round(loss.cpu().deatach().item(), 4)}.pt"
                                       )
                           )


def main(args):
    run = wandb.init(project='sxfst', config=args)
    save_path = os.path.join('weights', run.name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
               
    data = DataTensors(args.input, test=True, n_non_binders=3)
    #data = Data(args.input, test=True)
    a = len(data) // 4
    train_data, test_data = random_split(data, (len(data)-a, a))
    train_loader = DataLoader(train_data,
                              batch_size=args.batch_size,
                              shuffle=True,
                              num_workers=0,
                              )
    test_loader = DataLoader(test_data,
                             batch_size=args.batch_size,
                             shuffle=True,
                             num_workers=0,
                             )

    model = Model(esm=Esm(model=args.esm),
                  seqpool=SeqPool(conv_channels=35,
                                  num_conv_layers=args.num_conv_layers_pool,
                                  kernel_size=args.kernel_size_pool,
                                  stride=args.stride_pool,
                                  num_lstm_layers=args.num_lstm_layers_pool,
                                  lstm_hs=args.lstm_hs_pool,
                                  ),
                  fpnn=Fpnn(fp_size=2048,
                            emb_size=args.emb_size_fp,
                            n_layers=args.n_layers_fp,
                            ),
                  head=Head(emb_size=args.emb_size_head,
                            n_layers=args.n_layers_head,
                            )
                  )
    if args.cuda:
        model = model.cuda()

    wandb.watch(model)

    train(model,
          train_loader,
          proc=False,
          n_non_binders=1,
          lr=args.lr,
          epochs=args.epochs,
          save_path=save_path,
          cuda=args.cuda,
          )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-b', '--batch_size', default=2, type=int)
    parser.add_argument('-e', '--epochs', default=8, type=int)
    parser.add_argument('-l', '--lr', default=1e-6, type=float)
    parser.add_argument('--esm', default='esm1_t6_43M_UR50S')
    parser.add_argument('--cuda', action='store_true')
    # Head
    parser.add_argument('--emb_size_head', default=192, type=int)
    parser.add_argument('--n_layers_head', default=4, type=int)
    # Fpnn
    parser.add_argument('--emb_size_fp', default=128, type=int)
    parser.add_argument('--n_layers_fp', default=3, type=int)
    # SeqPool
    parser.add_argument('--num_conv_layers_pool', default=2, type=int)
    parser.add_argument('--num_lstm_layers_pool', default=2, type=int)
    parser.add_argument('--kernel_size_pool', default=9, type=int)
    parser.add_argument('--stride_pool', default=3, type=int)
    parser.add_argument('--lstm_hs_pool', default=32, type=int)
    args = parser.parse_args()
    assert args.input is not None
    assert os.path.isfile(args.input)
    main(args)
