#!/usr/bin/env python
import sys
import os
import re
import yaml
from pprint import pprint, pformat
import ast
from functools import lru_cache

import pandas as pd

from sxfst import utils

def find(path):
    ''' unix find
    '''
    return os.popen(f"find {path}").read().splitlines()

def grep(_list, regex):
    return list(filter(\
            lambda s : re.search(s, regex, re.IGNORECASE) is not None,
                _list))


class Config:
    def __init__(self,
                 find_results,
                 *args,
                 **kwargs,
                 ):
        self.find_results = find_results
        self.config_path = self.find_config(find_results)
        with open(self.config_path) as f:
            self.data = yaml.safe_load(f)
    def __getitem__(self, idx):
        return self.data[idx]
    def __repr__(self):
        return f"{pformat(self.data)}"
    def find_config(self, find_results):
        matches = [i for i in find_results if 'config' in i]
        assert len(matches) == 1 # future issue
        return matches[0]

class EchoFiles:
    def __init__(self,
                 find_results,
                 *args,
                 **kwargs,
                 ):
        self.find_results = find_results
        self.echo_files = [i for i in find_results if \
                           'echo' in i]
    def __len__(self):
        return len(self.find_results)
    def __getitem__(self,i):
        return self.echo_files[i]

class PlatereaderFiles:
    def __init__(self,
                 find_results,
                 *args,
                 **kwargs,
                 ):
        self.find_results = find_results

class Screen:
    ''' Container for screening data.
        Maps compounds to plate wells.
        Processes platereader data on a by compound basis.
    '''
    def __init__(self,
                 path,
                 *args,
                 **kwargs,
                 ):
        self.map_fs(path)
        self.map_wells()
    def __repr__(self):
        return f"{self.fs}"
    def map_fs(self, path):
        ''' Doesnt really return a map
            just finds echo, platereader and config files
            Sets attributes:
                self.fs 
                self.echo_files 
                self.platereader_files 
                self.config 
                self.picklists
        '''
        fs = find(path)
        self.fs = fs
        self.echo_files = EchoFiles(fs)
        self.platereader_files = PlatereaderFiles(fs)
        self.config = Config(fs)
        self.picklists = []
        for i in self.config['echo']['picklists']:
            for j in self.echo_files:
                # only checks if filenames match config.picklists
                if os.path.basename(i) == os.path.basename(j):
                    if i not in self.picklists:
                        self.picklists.append(j)
    def map_wells(self):
        ''' Maps compounds in picklists to platereader wells.
        '''
        if len(self.picklists) == 1:
            picklist_path = self.picklists[0]
            assert os.path.exists(picklist_path), \
                    f"path not found: {picklist_path}"
            picklist = proc_picklist(picklist_path)
            cpd_map = {}
            for cpd in set(picklist['Cpd']):
                o = {}
                chunk = picklist.loc[picklist['Cpd']==cpd, :]
                o['df'] = chunk
                echo_name_ = list(set(chunk['Destination Plate Name']))
                assert len(echo_name_) == 1
                echo_name = echo_name_[0]
                wells = list(chunk['DestWell'])
                m_ = map_echo_name_2path(echo_name, self.config)
                if m_ is not None:
                    m = m_[list(m_.keys())[0]]
                    if 'test' in m:
                        o['test'] = read_plate_csv(m['test'][0])[wells]
                    if 'control' in m:
                        o['ctrl'] = m['control']
                    if 'ctrl' in m:
                        o['ctrl'] = m['ctrl']
                    if 'echo_map' in m or 'echo_map' in m:
                        o['echo_map'] = m['echo_map']
                cpd_map[cpd] = o
            pprint(cpd_map)
        else:
            pass #print(self.picklists)

@lru_cache(128)
def read_plate_csv(path):
    return utils.PlateData(path)

def proc_picklist(path):
    picklist = pd.read_csv(path)
    if (unn:='Unnamed: 0') in picklist.columns:
        picklist = picklist.drop(unn, axis=1)
    # Cpd names are like : ['SXXXX']
    cpds_ = [ast.literal_eval(i) for i in picklist['Cpd']]
    cpds = [';'.join(i) for i in cpds_]
    picklist['Cpd'] = cpds 
    return picklist

@lru_cache(128)
def map_echo_name_2path(echo_name, config):
    platereader_files = config['platereader']
    for i,j in zip(platereader_files.keys(),
                   platereader_files.values()):
        item = j['echo_map']
        if item is not None:
            if echo_name == item[0]:
                return {i:j}

def proc(path):
    screen = Screen(path)
    #print(screen.config)
    #pprint(screen)

def main(args):
    for i in args:
        proc(i)

if __name__ == '__main__':
    main(sys.argv[1:])
