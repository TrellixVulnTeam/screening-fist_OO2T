#!/usr/bin/env python
import sys
import os
import re
import yaml
from pprint import pprint, pformat
import ast
from functools import lru_cache

import pandas as pd
from scipy.ndimage import gaussian_filter1d

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

class Cpd:
    def __init__(self,
                 **kwargs,
                 ):
        self.__dict__ = {**self.__dict__, **kwargs}
    def __repr__(self):
        return f"Cpd: {pformat(self.__dict__)}"

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
    def __getitem__(self, idx):
        assert hasattr(self, 'cpd_map')
        assert idx in self.cpd_map.keys()
        return self.cpd_map[idx]
    def __iter__(self):
        assert hasattr(self, 'cpd_map')
        for i in self.cpd_map:
            yield self[i]
    def __len__(self):
        assert hasattr(self, 'cpd_map')
        return len(self.cpd_map)
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
        self.root = path
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
                m = map_echo_name_2path(echo_name, 
                                        self.config,
                                        self.root,
                                       )
                if m is not None:
                   ctrl = read_plate_csv(m['ctrl'])[wells] if 'ctrl' in m.keys() else None
                   test = read_plate_csv(m['test'])[wells] if 'test' in m.keys() else None
                   echo_map = m['echo_map'] if 'echo_map' in m.keys() else None
                   cpd_map[cpd] = Cpd(name=cpd,
                                      test=test,
                                      ctrl=ctrl,
                                      echo_map=chunk,
                                      )
            self.cpd_map = cpd_map

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
def map_echo_name_2path(echo_name, config, root):
    add_root = lambda path : os.path.join(root, path)
    platereader_files = config['platereader']
    o = {}
    for i,j in zip(platereader_files.keys(),
                   platereader_files.values()):
        item = j['echo_map']
        if item is not None:
            if echo_name == item[0]:
                o['test'] = add_root(j['test'][0])
                o['control'] = add_root(j['control'][0])
                o['echo_map'] = add_root(j['echo_map'][0])
                return o

def norm_traces(test, 
                ctrl=None):
    '''
    '''
    test = test.subtract(test.loc[:,800],
                         axis=0)
    if ctrl is not None:
        ctrl = ctrl.subtract(ctrl.loc[:,800],
                             axis=0)
        norm = test - ctrl
        return norm
    else:
        return test

def smooth(df,
           sigma=3,
           ):
    cols = df.columns
    idx = df.index
    smth = pd.DataFrame(gaussian_filter1d(df, sigma),
                        index=idx,
                        columns=cols)
    return df

def is_anomaly(norm):
    return False

def diff(norm, baseline):
    return norm.subtract(baseline,
                         axis=0)

def response(_diff):
    a420 = _diff.loc[:,420]
    a390 = _diff.loc[:,390]
    return a420.abs().add(a390.abs())

def c2(v1, c1, v2):
    return (v1 * c1) / v2

def proc(cpd, # Cpd()
         sigma=3,
         ):
    name = cpd.name
    test = cpd.test
    ctrl = cpd.ctrl
    echo_map = cpd.echo_map
    norm = norm_traces(test, ctrl)
    smth = smooth(norm, sigma)
    if not is_anomaly(smth):
        change = diff(smth, smth[list(smth.keys())[0]])
        y = dict(response(change))
        #print(y)
        vol_well_map = dict(zip(echo_map['DestWell'],
                                echo_map['Transfer Volume /nl']))
        x = {i:c2(v1=vol_well_map[i] * 1e-9,    # nl to l
                  c1=10e-3,                     # 10 mM stock
                  v2=40e-6,                     # 40 ul well vol
                  ) / 1e-6                     # M to uM
              for i in vol_well_map}
        print(x)
        print(y)


def main(args):
    for path in args:
        screen = Screen(path)
        for cpd in screen:
            proc(cpd)

if __name__ == '__main__':
    main(sys.argv[1:])
