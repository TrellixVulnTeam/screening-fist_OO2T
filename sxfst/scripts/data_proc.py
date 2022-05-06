#!/usr/bin/env python
import sys
import os
import yaml
from pprint import pprint
import copy
from tqdm import tqdm                                                      
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
plt.style.use('dark_background')

from sxfst.utils import PlateData


class Config:
    def __init__(self, path=None, data=None):
        if path is not None:
            with open(path) as f:
                data = yaml.full_load(f)
        assert data is not None
        self.__dict__ = {**self.__dict__, **data}
    def __repr__(self):
        return yaml.dump(self.__dict__)
                            
    
def get_actual_vol(src_plate_name, 
                   dest_plate_name, 
                   exceptions,
                   well, 
                   vol):
    chunk = exceptions.loc[exceptions['src_plate'] == src_plate_name, :]
    chunk = chunk.loc[chunk['dest_plate'] == well, :]
    row = chunk.loc[chunk['Destination Well'] == well, :]
    if len(row) > 0:
        actual_vol = row['Actual Volume']
        return actual_vol
    if len(row) > 1:
        pass
    else:
        return vol
                    
def plot_set(*data, 
             labels=None,
             vols=None,
             titles=None,
             axs=[],
             save=None,
            ):                                                 
    assert len(data) > 0, 'empty set given'                                   
    if titles is None:
        titles = [None] * len(data)
    elif isinstance(titles, str):
        titles = [titles]*len(data)
    assert len(titles) == len(data), f'N titles({len(titles)}) != N plots ({len(data)})'
    n_plots = len(data) + len(axs)
    if n_plots < 8:                                                         
        fig, ax = plt.subplots(1, 
                               len(data) + len(axs),                                     
                               figsize=(4*n_plots,4))                                 
    else:                                                                      
        fig, ax = plt.subplots(n_plots//(nRows:=8),                         
                               nRows,                                          
                               figsize=(16,128*(len(data) / len(data))))     
                                                                               
    for plot_num, (wells, ax_, title) in enumerate(zip(data, ax.flatten(), titles)):                                  
        if vols is None:
            vols = range(len(wells))
        for row_, vol in zip(wells.index, vols):                                 
            ax_.plot(wells.loc[row_,:], 
                     c=plt.cm.cool(vol/2000),
                     label=vol)
                                                                               
        ax_.set_xlim(280,800)                                                  
        ax_.set_ylim(-0.02,0.25)                                                  
        ax_.set_title(title)                                                     
        ax_.set_xticks(wells.columns[::50])
        ax_.set_xlabel('Wavelength (nm)')                                      
        #ax_.axis('off')                                                        
        ax_.legend()
    if len(axs) > 0:
        for i,j in enumerate(axs, plot_num):
            ax_ = ax.flatten()[i]
            ax_ = j
                                                                               
    #plt.title(title)                                                           
    plt.tight_layout()                     
    if save is not None:
        assert isinstance(save, str)
        plt.savefig(save)
        plt.close()
    else:
        plt.show()
    
    
def norm_traces(test, 
                ctrl=None):
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
           axis=-1,
           ):
    cols = df.columns
    idx = df.index
    smth = pd.DataFrame(gaussian_filter1d(df, sigma, axis=axis),
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


def main(config_path):
    os.chdir(os.path.dirname(config_path))
    config = Config(os.path.basename(config_path))
    picklist = pd.read_csv(config.echo['picklist'][0], index_col=0)
    plate_names = dict(zip(picklist['Destination Plate Name'].unique(),
                           [f'plate_{i}' for i in range(1,16)]))

    picklist['dest_plate'] = [plate_names[i] for i in picklist['Destination Plate Name']]
    _exceptions = [i for i in  config.echo['transfers'] if 'Exceptions' in i]
    assert len(_exceptions) == 1
    _exceptions = _exceptions[0]

    exceptions = pd.read_csv(_exceptions)
    src_plate_names = dict(zip(sorted(exceptions['Source Plate Name'].unique()),
                               sorted(picklist['SrcPlate'].unique()),
                              ))
    dest_plate_names = dict(zip([f'Destination[{i}]' for i in range(2,17)],
                                 sorted(picklist['Destination Plate Name'].unique()),
                               ))

    exceptions['src_plate'] = [src_plate_names[i] for i in exceptions['Source Plate Name']]
    exceptions['dest_plate'] = [dest_plate_names[i] for i in exceptions['Destination Plate Name']]
    print(f'number of exceptions: {len(exceptions)}')
    a = [get_actual_vol(src_plate_name=i, 
                        dest_plate_name=j, 
                        well=j, 
                        vol=l,
                        exceptions=exceptions,
                        ) for i,j,k,l in zip(picklist['SrcPlate'],
                                                   picklist['Destination Plate Name'],
                                                   picklist['DestWell'],
                                                   picklist['Transfer Volume /nl'])]

    picklist['actual_vol'] = a

    plates = copy.deepcopy(config.platereader)
    for i in plates:
        plates[i]['picklist'] = picklist.loc[picklist['dest_plate'] == i, :]
        plates[i]['test']['data'] = PlateData(plates[i]['test']['path'])
        plates[i]['control']['data'] = PlateData(plates[i]['control']['path'])
        

    if not os.path.exists('tmp'):
        os.mkdir('tmp')

    sigma = 1/32

    for i in plates.keys():
        plate_i = plates[i]
        pck_i = plate_i['picklist']
        if 'control' in plate_i:
            #pprint(plate_i['control']['data'].__dict__)
            no_cpd = [i for i in plate_i['control']['data'].df.index \
                         if i not in pck_i['DestWell'].to_list()]
        else:
            no_cpd = [i for i in plate_i['test']['data'].df.index \
                         if i not in pck_i['DestWell'].to_list()]


        for j in tqdm(pck_i['Cpd'].unique()):
            pck_chunk = pck_i.loc[pck_i['Cpd'] == j, :]
            cpd_name = pck_chunk['Cpd'].unique()[0]
            
            ctrl_wells_raw = plate_i['control']['data'][pck_chunk['DestWell'].to_list()]
            test_wells_raw = plate_i['test']['data'][pck_chunk['DestWell'].to_list()]
            
            #test_wells_clipped = test_wells_raw.where(test_wells_raw < 1).dropna(axis=1)
            #ctrl_wells_clipped = ctrl_wells_raw.where(ctrl_wells_raw < 1).dropna(axis=1)
            
            test_wells_clipped = test_wells_raw.loc[:,300:800]
            ctrl_wells_clipped = ctrl_wells_raw.loc[:,300:800]
            
            test_wells_smooth = smooth(test_wells_clipped, sigma=sigma)
            ctrl_wells_smooth = smooth(ctrl_wells_clipped, sigma=sigma)
            
            test_wells_norm = norm_traces(test_wells_clipped)
            ctrl_wells_norm = norm_traces(ctrl_wells_clipped)
            
            #ctrl_wells = norm_traces(smooth(ctrl_wells_clipped, sigma=sigma))
            #test_wells = test_wells.where(test_wells < 1).dropna(axis=1)
            #ctrl_wells = ctrl_wells.where(ctrl_wells < 1).dropna(axis=1)
            
            ctrl_plate_row = set([i[0] for i in 
                        ctrl_wells_raw.index.to_list() + ctrl_wells_raw.index.to_list()])
            assert len(ctrl_plate_row) == 1
            ctrl_plate_row = list(ctrl_plate_row)[0]
            _blanks =  [i for i in no_cpd if ctrl_plate_row in i]
            blanks_raw = norm_traces(smooth(plate_i['control']['data'][_blanks], 
                                            sigma=sigma))
            blanks_clipped = blanks_raw.where(blanks_raw < 1).dropna(axis=1)
            blanks_mean =  blanks_clipped.mean(axis=0)
            
            test_plate_row = set([i[0] for i in 
                        test_wells_raw.index.to_list() + ctrl_wells_raw.index.to_list()])
            assert len(test_plate_row) == 1
            test_plate_row = list(test_plate_row)[0]
            _prots =  [i for i in no_cpd if test_plate_row in i]
            prots_raw = norm_traces(smooth(plate_i['test']['data'][_prots], 
                                           sigma=sigma))
            prots_clipped = prots_raw.where(prots_raw < 1).dropna(axis=1)
            prots_mean = prots_clipped.mean(axis=0)
            
            test_minus_ctrl_norm = test_wells_norm - ctrl_wells_norm
            test_traces_with_zero =  pd.concat([prots_mean, test_wells_norm.T], axis=1).T
            ctrl_traces_with_zero =  pd.concat([blanks_mean, ctrl_wells_norm.T], axis=1).T
            
            vols = [0] + list(pck_chunk['actual_vol']) # match up logs
            test_traces_with_zero.index = vols
            ctrl_traces_with_zero.index = vols
            plot_set(test_traces_with_zero, 
                     (a:=test_traces_with_zero - blanks_mean),
                     (b:=a - a.iloc[0,:]),
                     smooth(a.div(a.loc[:,390], axis=0), sigma=2) / 16,
                     vols=vols,
                     titles=['test_traces_with_zero', 
                             'minus mean blank',
                             'diff',
                             '?',
                            ],
                     save=os.path.join('tmp',f'{cpd_name}-specs.png')
                     )

if __name__ == '__main__':
    main(sys.argv[1])
