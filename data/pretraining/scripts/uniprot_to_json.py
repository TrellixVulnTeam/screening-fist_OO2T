#!/usr/bin/env python
import sys
import os
import json
from pprint import pformat
import re

FIELDS = ['ID',
          'AC',
          'DT',
          'DE',
          'GN',
          'OS',
          'OG',
          'OC',
          'OX',
          'OH',
          'RN',
          'RP',
          'RC',
          'RX',
          'RG',
          'RA',
          'RT',
          'RL',
          'CC',
          'DR',
          'PE',
          'KW',
          'FT']


def read(path):
    '''
    generator that iterates through .dat file, returning seperate entries
    '''
    with open(path,'r') as f:
        entry = ''
        for line in f:
            if line[0:2] != '//':
                entry = ''.join([entry,line])
            else:
                yield entry
                entry = ''

def read_stdin(stdin): 
    entry = ''
    for line in stdin:
        if line[0:2] != '//':
            entry = ''.join([entry,line])
        else:
            yield entry
            entry = ''

def get_seq(entry):
    start = re.search('SQ.+', entry).end() # check this
    return entry[start:].replace('\n','').replace(' ','')

def parse_entry(entry):
    data = {i:';'.join(re.findall(f'{i}\s+(.+)', entry)) for i in FIELDS}
    data['seq'] = get_seq(entry)
    return data

def main_stdin(stdin):
    for entry in read_stdin(stdin):
        data = parse_entry(entry)
        sys.stdout.write(json.dumps(data)+'\n')

def main(args):
    for arg in args:
        assert os.path.exists(arg)
        assert os.path.isfile(arg)
        for entry in read(arg):
            data = parse_entry(entry)
            sys.stdout.write(json.dumps(data)+'\n')



if __name__ == '__main__':
    if not sys.stdin.isatty():
        main_stdin(sys.stdin)
    else:
        main(sys.argv[1:])
