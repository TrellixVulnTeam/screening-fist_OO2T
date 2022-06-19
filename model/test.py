#!/usr/bin/env python
import os
import argparse
from tqdm import tqdm
import json
import pickle as pkl

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from torch import Tensor, cat, mean
from einops import rearrange, repeat
from sklearn.metrics import confusion_matrix, precision_recall_curve, roc_curve, det_curve, average_precision_score

from data import Data, DataTensors, DataEmbeddings2
from model import Model, Model2, Esm, Head, Fpnn, SeqPool


def test(model,
         data_loader,
         cuda=False,
         n_non_binders=1,
         output=None,
         ):
    loss_fn = nn.BCELoss()#weight=weight)
    losses, ys, yhs = [], [], []
    with torch.no_grad():
        for seq_, fingerprints_, hit_ in tqdm(data_loader):
            if cuda:
                seq_, fingerprints_, hit_ = seq_.cuda(), fingerprints_.cuda(), hit_.cuda()
            #seq = rearrange(seq_, 'b1 b2 l -> (b1 b2) l')
            seq = rearrange(seq_, 'b1 b2 l d -> (b1 b2) l d') # embeddings
            fingerprints = rearrange(fingerprints_, 'b1 b2 l -> (b1 b2) l')
            hit = rearrange(hit_, 'b1 b2 -> (b1 b2)')
            yh = model(seq, fingerprints)
            if len(hit.shape) == 1:
                y = rearrange(hit.float(), '(b c) -> b c', c=1)
            else:
                y = hit.float()
            loss = loss_fn(yh, y)
            losses.append(loss.reshape(-1).detach())
            ys.append(y.reshape(-1).detach())
            yhs.append(yh.reshape(-1).detach())
    losses = cat(losses)
    mean_loss = mean(losses).cpu().numpy()
    ys = cat(ys)
    yhs = cat(yhs)
    yhs_ = (yhs > 0.5).cpu().numpy()
    ys = ys.cpu().numpy()
    yhs = yhs.cpu().numpy()


    cfz = confusion_matrix(ys, yhs_)
    (tp, fn), (fp, tn) = cfz
    precision, recall, pr_thresholds = precision_recall_curve(ys, yhs)
    fpr, tpr, roc_thresholds = roc_curve(ys, yhs)
    fprd, fnrd, thresholdsd = det_curve(ys, yhs)
    avpr = average_precision_score(ys, yhs_)
    d = {'mean_loss':mean_loss.tolist(),
         'average_precision_score':avpr.tolist(),
         'confusion_matrix':{'true_pos':tp.tolist(),
                             'true_neg':tn.tolist(),
                             'false_pos':fp.tolist(),
                             'false_neg':fn.tolist(),
                             },
         'roc':{'false_pos':fpr.tolist(),
                'true_pos':tpr.tolist(),
                'roc_thresholds':roc_thresholds.tolist(),
                },
         'precision_recall_curve':{'precision':precision.tolist(),
                                   'recall':recall.tolist(),
                                   'pr_thresholds':pr_thresholds.tolist(),
                                   },
         'det':{'fpr':fprd.tolist(),
                'fnr':fnrd.tolist(),
                'thresholds':thresholdsd.tolist()},
         }
    if not os.path.exists(output):
        os.makedirs(output)
    with open(os.path.join(output, 'val.json'), 'w') as f:
        json.dump(d, f)

def main(args):
    data = DataEmbeddings2(args.input, 
                           embeddings_dir='embeddings', 
                           max_seq_len=800,
                           n_non_binders=1,
                           cuda=args.cuda,
                           )
    data_loader = DataLoader(data,
                             batch_size=args.batch_size,
                             shuffle=True,
                             num_workers=8,
                             )

    with open(args.pkl, 'rb') as f:
        model = pkl.load(f)

    if args.cuda:
        model = model.cuda()


    test(model=model,
         data_loader=data_loader,
         n_non_binders=1,
         cuda=args.cuda,
         output=args.output,
         )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-p', '--pkl')
    parser.add_argument('-o', '--output')
    parser.add_argument('-b', '--batch_size', default=128, type=int)
    parser.add_argument('--cuda', action='store_true')
    args = parser.parse_args()
    assert args.pkl is not None
    assert os.path.isfile(args.pkl)
    assert args.input is not None
    assert os.path.isfile(args.input)
    main(args)

