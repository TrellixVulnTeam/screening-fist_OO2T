import argparse
from tqdm import tqdm
import swalign

def parse(path):
    # assumes one sequence
    with open(path,'r') as f:
        data=f.read().splitlines()
    header = data[0]
    seq = ''.join(data[1:])
    return header, seq

def main(args):

    match=2
    mismatch=-1
    scoring=swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    
    wt_head, wt_seq = parse(args.ref)
    wt_seq = wt_seq.upper()

    with open(args.out,'w') as f:
        for i in tqdm(args.sequences):
            head, seq = parse(i)
            aln = sw.align(wt_seq, seq)
            f.write(f'{head} and {wt_head} \n')
            aln.dump(wrap = 60, out=f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--out', default='out')
    parser.add_argument('-r','--ref')
    parser.add_argument('-s','--sequences', nargs='+')
    args = parser.parse_args()
    main(args)
