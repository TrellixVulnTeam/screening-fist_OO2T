#!/usr/bin/env python
import sys
from functools import lru_cache
from esm import Alphabet
from rdkit import Chem



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

abc = Alphabet(TOKS)

def filter_seq(seq):
    try: 
        abc.encode(seq)
        return True
    except Exception as e:
        #print(e)
        return False

@lru_cache(128)
def filter_smiles(smiles):
    try:
        Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))
        return True
    except Exception as e:
        #print(e)
        return False

def main(args):
    for arg in args:
        with open(arg) as f:
            for i, line in enumerate(f):
                if i == 0:
                    sys.stdout.write(line)
                seq, smiles, hit = line.split(',')
                if filter_seq(seq) and filter_smiles:
                    sys.stdout.write(line)
    

def mainin(stdin):
    header = stdin.readline()
    sys.stdout.write(header)
    for line in stdin:
        seq, smiles, hit = line.split(',')
        if filter_seq(seq) and filter_smiles:
            sys.stdout.write(line)

if __name__ == '__main__':
    if not sys.stdin.isatty():
        mainin(sys.stdin)
    else:
        main(sys.argv[1:])
