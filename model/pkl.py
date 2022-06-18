#!/usr/bin/env python
import os
import argparse
import pickle as pkl
import torch
from model import Model2, SeqPool, Fpnn, Head

def main(args):
    model = Model2(seqpool=SeqPool(conv_channels=35,
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
                            layer={True:'transformer', False:'linear'}[args.transformer],
                            )
                  )
    model.load_state_dict(torch.load(args.weights))

    with open(os.path.join(os.path.dirname(args.weights), 
                           os.path.basename(args.weights).split('.')[0] +'.pkl'),
              'wb') as f:
        pkl.dump(model, f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--weights', help='path to weights')
    parser.add_argument('-b', '--batch_size', default=2, type=int)
    parser.add_argument('-e', '--epochs', default=8, type=int)
    parser.add_argument('-l', '--lr', default=1e-6, type=float)
    parser.add_argument('--esm', default='esm1_t6_43M_UR50S')
    parser.add_argument('--cuda', action='store_true')
    parser.add_argument('--wandb', action='store_true')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--transformer', action='store_true')
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
    parser.add_argument('--load')
    args = parser.parse_args()
    assert args.weights is not None
    assert os.path.isfile(args.weights)
    main(args)
