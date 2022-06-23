#!/usr/bin/env python
import sys
import argparse
import os
import pickle as pkl
import random
import heapq
from tqdm import tqdm

import pandas as pd
from rdkit import Chem
import torch
from torch import cat, Tensor, zeros
from einops import rearrange, repeat

import esm
from esm import Alphabet
from model import PredModel, fp
from data import TOKS

AAS = list('ACDEFGHIKLMNPQRSTVWY')

MXN_SITES = [47, 49, 51, 75, 78, 88, 94, 138, 142, 175, 178, 184, 188, 
             205, 226, 252, 255, 260, 263, 290, 295, 328, 330, 350, 353]

MESOTRIONE = 'CS(=O)(=O)C1=CC(=C(C=C1)C(=O)C2C(=O)CCCC2=O)[N+](=O)[O-]'

BM3_A82F = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQALKFVRDFFGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK'

def mutate(seq,
           pos:int,
           new:str):
    ''' mutate string at pos to new
    '''
    seq = list(seq)
    seq[pos] = new
    return ''.join(seq)

def random_mutate(seq:str,
                  vocab=AAS,
                  pos_weights=None,
                  vocab_weights=None):
    ''' mutate string at random position to random characters
        from vocab (iterable of chars),
        seq : str
        vocab : iterable returning strings
        pos_weights : iterable of floats, maps to seq
            probability weights for position selection
        vocab_weights : iterable of floats, maps to vocab
            probability weights for substitution selection
    '''
    mxn_site = random.choices(range(len(seq)), weights=pos_weights, k=1)[0]
    mxn_site = random.choices(MXN_SITES, weights=pos_weights, k=1)[0]
    new = random.choices(vocab, weights=vocab_weights, k=1)[0]
    return mutate(seq, mxn_site, new)

def crossover(a:str,
              b:str):
    ''' randomly splice two strings
        returns string
    '''
    cut = random.randint(1,min(len(a),len(b))-1)
    return random.choice([a[:cut] + b[cut:], b[:cut] + a[cut:]])

def mutate_string(template, target_dict):
    s_ = list(template)
    for i,j in zip(target_dict.keys(), target_dict.values()):
        s_[i] = j
    return ''.join(s_)

def main(args):
    max_seq_len = 800
    if args.smiles is None:
        smiles = MESOTRIONE
    else:
        smiles = args.smiles
    if args.sequence is not None:
        sequence_ = args.sequence
    else:
        sequence_ = BM3_A82F

    model = PredModel(args.pkl).train()
    model.eval()

    alphabet = Alphabet(TOKS)
    pad = lambda x, l : cat([x, zeros(l-len(x))])
    encode_seq = lambda seq : pad(Tensor(alphabet.encode(seq)), max_seq_len).unsqueeze(0)

    ## ga loop
    fp_tensor = fp(smiles)
    if args.cuda:
        model = model.cuda()
        fp_tensor = fp_tensor.cuda()

    fp_tensor_batch = repeat(fp_tensor, 'b l -> (b b2) l', b2=args.pop_size)
    #seqx_batch = repeat(seqx, 'b l -> (b b2) l', b2=args.pop_size)
    o_ = pd.DataFrame([], columns=['generation', 'score', 'sequence'])
    o_.to_csv(sys.stdout, index=False)
    
    pop = [random_mutate(sequence_) for _ in range(args.pop_size)]
    bestest = {}
    with tqdm(range(args.epochs)) as bar:
        for i, _ in enumerate(bar):
            seqx_batch = cat([encode_seq(i) for i in pop], dim=0)
            if args.cuda:
                seqx_batch = seqx_batch.cuda()
            yh = model(seqx_batch, fp_tensor_batch).detach().cpu()
            pop_scores = dict(zip(pop, yh))
            best = heapq.nlargest(args.pop_size // 4, 
                                  {**pop_scores, **bestest},
                                  key=lambda i : {**pop_scores, **bestest}[i])
            bestest = {i:{**pop_scores, **bestest}[i] for i in best}

            bar.set_postfix({'max':max(pop_scores.values()).item(),
                             'mean':yh.mean().item(),
                             'min':min(pop_scores.values()).item(),
                             })
            bestest_ = heapq.nlargest(args.pop_size, 
                                      bestest,
                                      key=lambda i : bestest[i])
            bestest = {i:bestest[i] for i in bestest_}

            pop = [random_mutate(\
                    crossover(*random.choices(list(bestest.keys()), k=2))) \
                        for _ in range(args.pop_size)]
        
            o = pd.DataFrame({'generation':[i]*len(pop_scores),
                               'scores':[i.item() for i in pop_scores.values()],
                               'sequence':pop_scores.keys()},
                             )

            #o = o.loc[:,['generation', 'score', 'sequence']]
            o.to_csv(sys.stdout, index=False, header=False)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--smiles', help='smiles code to optimize sequence to')
    parser.add_argument('-p', '--pkl', help='model pickle file')
    parser.add_argument('-S', '--sequence')
    parser.add_argument('-n', '--pop_size', default=32, type=int)
    parser.add_argument('-e', '--epochs', default=32, type=int)
    parser.add_argument('--cuda', action='store_true')
    args = parser.parse_args()
    assert args.pkl is not None
    assert os.path.isfile(args.pkl)
    main(args)

