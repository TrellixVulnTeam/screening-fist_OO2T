#!/usr/bin/env python
import sys
import os
from sxfst.utils import PlateData

def get_id(path):
    pass

def main(args):
    for arg in args:
        if os.path.exists(args):
            if os.path.isdir(arg):
                find = os.popen(f'find {arg}').read().split('\n')
                for _file in find:
                    print(get_id(_file))
            elif os.path.isfile(arg):
                print(get_id(_file))
        else:
            raise Warning(f'OS Error: {arg} not found')

if __name__ == '__main__':
    main(sys.argv[1:])
