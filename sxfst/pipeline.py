import sys
import os
import re
import json
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

        todo: 
            - echo exceptions report
            - get baselines
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
        if idx in self.cpd_map.keys():
            return self.cpd_map[idx]
        #elif idx in self.plates
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
        # else concat
        else:
            raise Warning('time to make a concat picklists thing')
        cpd_map = {}                       # becomes self.cpd_map, keys: cpd names
                                           # values: Cpd data object
        for cpd in set(picklist['Cpd']):
            o = {}
            chunk = picklist.loc[picklist['Cpd']==cpd, :]
            o['df'] = chunk
            dest_plate_ = list(set(chunk['Destination Plate Name']))
            assert len(dest_plate_) == 1
            dest_plate = dest_plate_[0]
            wells = list(chunk['DestWell'])
            # returns dict with keys:
            # 'test'     : path to test plate
            # 'control'  : path to ctrl plate
            # 'echo_map' :
            m = map_echo_name_2path(dest_plate, 
                                    self.config,
                                    self.root,
                                   )
            if m is not None:
               if 'ctrl' in m.keys():
                   ctrl_plate = read_plate_csv(m['ctrl'])
               elif 'control' in m.keys():
                   ctrl_plate = read_plate_csv(m['control'])
               else:
                   ctrl_plate = None

               # lru_cache fn, return: utils.Plate objects
               test_plate = read_plate_csv(m['test']) if 'test' in m.keys() else None

               # traces for ctrl and test wells
               ctrl = ctrl_plate[wells] if ctrl_plate is not None else None
               test = test_plate[wells] if test_plate is not None else None

               # get no compound wells from same row for test and ctrl plates
               # same row because of multichannel pipetting direction
               _picklist = picklist.loc[picklist['Destination Plate Name'] == dest_plate,:]

               if ctrl is not None: 
                   ctrl_no_cpd_wells = list(set(ctrl_plate.df.index).difference(\
                                                set(_picklist['DestWell'])))
                   baselines = {i:list(filter(lambda s : s[0] == i[0],
                                              ctrl_no_cpd_wells))
                                        for i,j in zip(ctrl.index, 
                                                       ctrl_no_cpd_wells)}
                   ctrl_baseline = {i:ctrl_plate[baselines[i]] for i in ctrl.index}
               else:
                    ctrl_baseline = None

               if test is not None:
                   test_no_cpd_wells = list(set(test_plate.df.index).difference(\
                                                set(_picklist['DestWell'])))
                   baselines = {i:list(filter(lambda s : s[0] == i[0],
                                              test_no_cpd_wells))
                                        for i,j in zip(test.index, 
                                                       test_no_cpd_wells)}
                   test_baseline = {i:test_plate[baselines[i]] for i in test.index}
               # v 
               echo_map = m['echo_map'] if 'echo_map' in m.keys() else None
               smiles = self.lib.loc[cpd, 'SMILES'] # use case-specific
               cpd_map[cpd] = Cpd(name=cpd,
                                  test=test,
                                  test_path=m['test'] if 'test' in m.keys() else None,
                                  ctrl=ctrl,
                                  ctrl_path=m['ctrl'] if 'ctrl' in m.keys() else None,
                                  ctrl_baseline=ctrl_baseline,
                                  test_baseline=test_baseline,
                                  echo_map=chunk,
                                  smiles=smiles,
                                  )
        self.cpd_map = cpd_map

    def map_baseline_wells(self):
        ''' find wells with no compounds, for baselines
        '''
        pass


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
def map_echo_name_2path(dest_plate, config, root):
    ''' map echo plate name to platereader data
        from config
    '''
    add_root = lambda path : os.path.join(root, path)
    platereader_files = config['platereader']
    o = {}
    for i,j in zip(platereader_files.keys(),
                   platereader_files.values()):
        item = j['dest_plate']
        if item is not None:
            if dest_plate == item[0]:
                o['test'] = add_root(j['test'][0])
                o['control'] = add_root(j['control'][0])
                o['dest_plate'] = j['dest_plate'][0]
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

def diff(norm, baseline, axis=0):
    ''' norm - baseline (axis=0)
    '''
    return norm.subtract(baseline.reset_index(drop=True),
                         axis=axis,
                         ).fillna(0)

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
        x_ = np.nan_to_num(np.array(list(x.values())), nan=0)
    if isinstance(y,dict):
        y_ = np.nan_to_num(np.array(list(y.values())), nan=0)
    try:
        (km, vmax), covariance = curve_fit(michaelis_menten, x_, y_, 
                bounds=((0, 0),
                        (1e4, 0.1),
                        ))
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
    ctrl_baseline = cpd.ctrl_baseline
    test_baseline = cpd.test_baseline
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
    o['test'] = {'path':cpd.test_path,
                 'data':test.to_json() if test is not None else test,
                 }
    o['ctrl'] = {'path':cpd.ctrl_path,
                 'data':ctrl.to_json() if ctrl is not None else ctrl,
                 }
    o['echo_map'] = echo_map.to_json() if echo_map is not None else echo_map
    o['x'] = x
    o['smiles'] = smiles
    norm = norm_traces(test, ctrl)
    smth = smooth(norm, sigma)

    # get protein baseline, might break
    _test_baseline = pd.concat(test_baseline.values()).drop_duplicates()
    _normbt = norm_traces(_test_baseline) # protein only wells
    _smthbt = smooth(_normbt, sigma) # smooth protein only traces
    prot_baseline = _smthbt.mean(axis=0)

    if not is_anomaly(smth):
        change = diff(smth, 
                      prot_baseline)
        y = dict(response(change))
        mm = fit_michaelis_menten(x,y) # {'km':km, 'vmax':vmax, 'rsq':rsq}
        o['y'] = y
        o['michaelis_menten'] = mm
        if plot:
            utils.plot_report(traces=smth,
                              diff=change,
                              name=name,
                              x=x,
                              y=y,
                              mm=mm,
                              save_path=os.path.join('img', name),
                              smiles=smiles,
                              )
    return o 


def main(args):
    for path in args:
        screen = Screen(path)
        #from multiprocessing.pool import ThreadPool
        #def helper(cpd):
        #    return proc(cpd, plot=True)
        #with ThreadPool(16) as pool: # pyplot 20 open warning
        #    pool.map(helper, 
        #             screen)
        #    pool.join()
        #for cpd in screen:
        for cpd in tqdm(screen):
            data = proc(cpd, plot=True)
        #pprint(data)

if __name__ == '__main__':
    main(sys.argv[1:])
