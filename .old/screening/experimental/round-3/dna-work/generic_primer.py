import re
from tqdm import tqdm 
import numpy as np
import pandas as pd 
from scipy import optimize

import mxn
from mxn.agilent_tm import agilent_tm
from bm3 import bm3, orf

import primer3

COMPLIMENT = {'A':'T','T':'A','C':'G','G':'C'}

reverse_comp = lambda s : ''.join([COMPLIMENT[i.upper()] for i in s[::-1]])


def agilent_tm(s):
    mv_conc, dv_conc, dntp_conc, dna_conc, dmsoPerc = [25, # mv - 25-100
            0.5, # dv - 0.5 - 5
            0.8, # dntp
            1, # dna
            8] # dmso 0-10%
    tm = primer3.bindings.calcTm(s, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc,
            tm_method = 'santalucia',
            salt_corrections_method = 'santalucia')
    return tm - 0.75 * dmsoPerc

def generic_primer(seq, comp='rev', tm=78):
    # seq = sense, region to make primer in
    if comp == 'rev':
        seq = reverse_comp(seq)

    homoTm = lambda s : primer3.calcHomodimerTm(s,mv_conc= 25, dv_conc = 0.5)
    endscore = lambda s : sum([1 for i in s[1] + s[-2] if i == 'C' or i == 'G']\
                            + [2 for i in s[0] + s[-1] if i == 'C' or i == 'G']) # max 6
    scorefn = lambda s : abs(agilent_tm(s) - tm) \
                        - (endscore(s) / 6) \
                        + max([(homoTm(s) / 6),0])
    select = lambda n1, n2 : seq[round(n1):round(n1)+round(n2)]
    objective = lambda n1, n2 : scorefn(select(n1, n2))
    helper = lambda array : objective(array[0], array[1])
    
    results = optimize.dual_annealing(helper, 
                            bounds = ((0,len(seq) - 60),
                                    (10,60)),
                            initial_temp = 1e5)
    primer = select(results.x[0], results.x[1])
    score = scorefn(primer)
    tm = agilent_tm(primer)
    return {'primer':primer,
            'tm':tm,
            'end_score': endscore(primer),
            'homotm':homoTm(primer),
            'length':len(primer)}

def main():
    primer = generic_primer(orf[-100:])





if __name__ == '__main__':
    main()
