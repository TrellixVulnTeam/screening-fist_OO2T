import sys
import os
import re
import yaml
from pprint import pprint, pformat
import ast
from functools import lru_cache
from tqdm import tqdm

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

from sxfst import utils

def find(path):
    ''' unix find
    '''
    return os.popen(f"find {path}").read().splitlines()

def grep(_list, regex):
    ''' unix grep, not used(?)
    '''
    return list(filter(\
            lambda s : re.search(s, regex, re.IGNORECASE) is not None,
                _list))


class Config:
    ''' finds config.yml in list of paths,
        parses. behaves like dict
    '''
    def __init__(self,
                 find_results,
                 *args,
                 **kwargs,
                 ):
        self.find_results = find_results
        self.config_path = self.find_config(find_results)
        with open(self.config_path) as f:
            self.data = yaml.safe_load(f)
    @property
    def keys(self):
        return list(self.data.keys())
    def __getitem__(self, idx):
        return self.data[idx]
    def __repr__(self):
        return f"{pformat(self.data)}"
    def find_config(self, find_results):
        matches = [i for i in find_results if 'config' in i]
        assert len(matches) == 1 # future issue
        return matches[0]

class EchoFiles:
    ''' Container for echo_files, 
        behaves like list.
    '''
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

class Data:
    ''' data container 
    '''
    def __init__(self,
                 **kwargs,
                 ):
        self.__dict__ = {**self.__dict__, **kwargs}

class Cpd(Data):
    ''' data container for compounds
    '''
    def __init__(self,
                 **kwargs,
                 ):
        super().__init__(**kwargs)
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
        self.root = path
        self.map_fs(path)
        self.lib = self.get_lib()
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
    def get_lib(self):
        if hasattr(self, 'config'):
            lib_path_ = self.config['lib']
            assert len(lib_path_) == 1, 'two libs, time to fix this'
            lib_path = os.path.join(self.root, lib_path_[0])
            assert os.path.exists(lib_path)
            lib = pd.read_csv(lib_path)
            unn ='Unnamed: 0'
            if unn in lib.columns:
                lib = lib.drop(unn, axis=1)
            #if (unn:='Unnamed: 0') in lib.columns:
            #    lib = lib.drop(unn, axis=1)
            lib.index = lib['CatalogNumber'] # use case-specific
            return lib
        else:
            raise Warning('no lib')
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
                   smiles = self.lib.loc[cpd, 'SMILES'] # use case-specific
                   cpd_map[cpd] = Cpd(name=cpd,
                                      test=test,
                                      ctrl=ctrl,
                                      echo_map=chunk,
                                      smiles=smiles,
                                      )
            self.cpd_map = cpd_map


@lru_cache(128)
def read_plate_csv(path):
    return utils.PlateData(path)

def proc_picklist(path):
    ''' reads my weird picklists, return dataframe
    '''
    picklist = pd.read_csv(path)
    unn = 'Unnamed: 0'
    if unn in picklist.columns:
        picklist = picklist.drop(unn, axis=1)
    #if (unn:='Unnamed: 0') in picklist.columns:
    #    picklist = picklist.drop(unn, axis=1)
    # Cpd names are like : ['SXXXX']
    cpds_ = [ast.literal_eval(i) for i in picklist['Cpd']]
    cpds = [';'.join(i) for i in cpds_]
    picklist['Cpd'] = cpds 
    return picklist

@lru_cache(128)
def map_echo_name_2path(echo_name, config, root):
    ''' map echo plate name to platereader data
        from config
    '''
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
    ''' translate y so that 800x = 0
        subtract control traces from test 
        if supplied
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
    ''' gaussian smoothing, returns df
    '''
    cols = df.columns
    idx = df.index
    smth = pd.DataFrame(gaussian_filter1d(df, sigma),
                        index=idx,
                        columns=cols)
    return df

def is_anomaly(norm):
    ''' detect:
            - 
    '''
    return False

def diff(norm, baseline):
    ''' norm - baseline (axis=0)
    '''
    return norm.subtract(baseline,
                         axis=0)

def response(_diff):
    ''' changes in absorbance at key wavelengths
        sum
    '''
    a420 = _diff.loc[:,420]
    a390 = _diff.loc[:,390]
    return a420.abs().add(a390.abs())

def michaelis_menten(x, km, vmax):
    return ((x * vmax) / (km + x)) 


def r_squared(yi,yj):
    residuals = yi - yj
    sum_sq_residual = sum(residuals ** 2)
    sum_sq_total = sum((yi - yi.mean()) ** 2) # check this!!!
    return 1 - (sum_sq_residual / sum_sq_total)

def fit_michaelis_menten(x,y):
    if isinstance(x,dict):
        x_ = np.array(list(x.values()))
    if isinstance(y,dict):
        y_ = np.array(list(y.values()))
    try:
        (km, vmax), covariance = curve_fit(michaelis_menten, x_, y_, 
                bounds=((0, 0),(1e2,0.2)))
    except RuntimeError:
        km, vmax = np.inf, np.inf
    yh = michaelis_menten(x_, km, vmax)
    rsq = r_squared(y_, yh)
    return {'km':km, 'vmax':vmax, 'rsq':rsq}

def c2(v1, c1, v2):
    return (v1 * c1) / v2

def proc(cpd, # Cpd()
         baseline=None,
         sigma=3,
         plot=False,
         ):
    o = {} # output
    name = cpd.name
    test = cpd.test
    ctrl = cpd.ctrl
    smiles = cpd.smiles
    echo_map = cpd.echo_map
    vol_well_map = dict(zip(echo_map['DestWell'],
                            echo_map['Transfer Volume /nl']))
    x = {i:round(c2(v1=vol_well_map[i] * 1e-9, # nl to l
                 c1=10e-3,                     # 10 mM stock
                 v2=40e-6,                     # 40 ul well vol
                 ) / 1e-6,                     # M to uM
                 3)                            # round to 3 places
              for i in vol_well_map}
    o['name'] = name
    o['test'] = test
    o['ctrl'] = ctrl
    o['echo_map'] = echo_map
    o['x'] = x
    o['smiles'] = smiles
    norm = norm_traces(test, ctrl)
    smth = smooth(norm, sigma)
    
    if not is_anomaly(smth):
        change = diff(smth, 
                      smth[list(smth.keys())[0]])
        y = dict(response(change))
        mm = fit_michaelis_menten(x,y) # {'km':km, 'vmax':vmax, 'rsq':rsq}
        o['y'] = y
        o['michaelis_menten'] = mm
        if plot:
            utils.plot_report(traces=smth,
                              #diff=change,
                              name=name,
                              x=x,
                              y=y,
                              mm=mm,
                              save_path=os.path.join('img', name),
                              smiles=smiles,
                              )
                              #**o)
    #if plot:
    #    utils.plotTraces(smth,
    #                     save_path=name,
    #                     cpd=name,
    #                     vols=x.values(),  
    #                     save=True,
    #                     )
    #    o['plot'] = save_path

    return o 


def main(args):
    for path in args:
        screen = Screen(path)
        #from multiprocessing.pool import ThreadPool
        #def helper(cpd):
        #    return proc(cpd, plot=True)
        #with ThreadPool(32) as pool:
        #    pool.map(helper, 
        #             screen)
        #    pool.join()
        #for cpd in screen:
        for cpd in tqdm(screen):
            proc(cpd, plot=True)

if __name__ == '__main__':
    main(sys.argv[1:])
