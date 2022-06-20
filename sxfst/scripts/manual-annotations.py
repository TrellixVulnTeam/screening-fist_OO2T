#!/usr/bin/env python
import sys
import os
import re
import json

with open('smiles_chemid.json') as f:
    CHEMMAP = json.load(f)

MAP = {'01.0/': 'BM3 Heme WT', 
       '02.0/': 'BM3 Heme A82F/F87V', 
       '03.0/': 'BM3 Heme A82F', 
       '04.0/': 'BM3 Heme 1YQP', 
       '05.0/': 'BM3 Heme 1YQO'}

WT = 'MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK'

def mutate(seq, pos, replace):
    seq_ = list(seq)
    seq_[pos] = replace
    return ''.join(seq_)
    
MUTANTS = {'WT':WT,
           'A82F/F87V':mutate(mutate(WT,82,'F'), 87, 'V'),
           'A82F':mutate(WT,82,'F'),
           '1YQP':mutate(WT,268, 'N'),
           '1YQO':mutate(WT,268, 'A'),
           }

def search(p, s):
    srch = re.search(p,s)
    if srch is not None:
        return srch.group()

def proc(path):
    name = os.path.basename(path).split('.')[0]
    with open(path) as f:
        data = f.read().split('\n')
    o = []
    for i in data:
        if i != '':
            experiment = search('([0-9]+\.0/)',i)
            chemid = search('S[0-9]{4}',i)
            protein_name = MAP[experiment] 
            o.append({'experiment':experiment,
                      'protein':protein_name,
                      'seq':MUTANTS[protein_name.split()[-1]],
                      'chemid':chemid,
                      'smiles':CHEMMAP[chemid] ,
                      'annotation':name, # duplicates
                      'hit':True if name == 'hits' else False,
                      })
    return o


def main(args):
    fields = ['experiment', 
              'protein', 
              'seq', 
              'chemid', 
              'smiles', 
              #'annotation',
              'hit',
              ]
    tab = '\t'
    sys.stdout.write(f"{tab.join(fields)}\n")
    o = []
    for arg in args:
        if os.path.isdir(arg):
            for i in os.listdir(arg):
                path = os.path.join(arg, i)
                if os.path.isfile(path) and 'txt' in path:
                    data = proc(path)
                    for i in data:
                        o.append(\
                    f"{tab.join([str(i[j]) for j in fields])}\n")
        elif os.path.isfile(arg):
            data = proc(arg)
            for i in data:
                o.append(\
                    f"{tab.join([str(i[j]) for j in fields])}\n")
    for i in set(o):
        sys.stdout.write(i)


if __name__ == '__main__':
    main(sys.argv[1:])
