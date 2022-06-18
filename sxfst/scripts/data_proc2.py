#!/usr/bin/env python
import sys
import os
import re
import yaml
from pprint import pprint
import copy
from tqdm import tqdm                                                      
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
plt.style.use('dark_background')

from sxfst.utils import PlateData, find, grep
from sxfst.data import smooth, is_anomaly, diff, c2

'''
Data processing script,
reads experminet run config, picks out experiments and plots,
also outputs combined traces to a single csv in stdout for analysis
'''

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
                    

def get_dispense_map(picklist_path, exceptions_paths=[]):
    picklist = pd.read_csv(picklist_path, index_col=0)
    plate_names = dict(zip(picklist['Destination Plate Name'].unique(),
                           [f'plate_{i}' for i in range(1,16)]))

    dispense_map = picklist.copy()
    dispense_map['dest_plate'] = [plate_names[i] for i in dispense_map['Destination Plate Name']]
    _exceptions = [i for i in  exceptions_paths if 'Exceptions' in i]
    assert len(_exceptions) == 1
    _exceptions = _exceptions[0]

    exceptions = pd.read_csv(_exceptions)
    src_plate_names = dict(zip(sorted(exceptions['Source Plate Name'].unique()),
                               sorted(dispense_map['SrcPlate'].unique()),
                              ))
    dest_plate_names = dict(zip([f'Destination[{i}]' for i in range(2,17)],
                                 sorted(dispense_map['Destination Plate Name'].unique()),
                               ))

    exceptions['src_plate'] = [src_plate_names[i] for i in exceptions['Source Plate Name']]
    exceptions['dest_plate'] = [dest_plate_names[i] if i in dest_plate_names \
            else re.search('Destination\[[0-9]+\]', i).group()
            for i in exceptions['Destination Plate Name']]
    a = []
    for i,j,k,l in zip(dispense_map['SrcPlate'],
                       dispense_map['Destination Plate Name'],
                       dispense_map['DestWell'],
                       dispense_map['Transfer Volume /nl']):
        a.append(get_actual_vol(src_plate_name=i, 
                            dest_plate_name=j, 
                            well=j, 
                            vol=l,
                            exceptions=exceptions,
                            ))

    dispense_map['actual_vol'] = a
    return dispense_map

def pad_trace(trace, 
              minn=220, 
              maxx=800,
              ):
    assert isinstance(trace, pd.Series)
    trace_wavelengths = trace.index
    padding_wavelengths = sorted(set(range(minn,maxx+1)).\
                                     difference(set(trace_wavelengths)))
    padding = pd.Series([np.nan] * len(padding_wavelengths),
                        index=padding_wavelengths,
                        dtype='float64')
    padded_trace = trace.append(padding).sort_index()
    return padded_trace

def proc_plate(plate,
               output_path,
               dispense_map_chunk,
               append=True,
               **kwargs
               ):
    bag_of_df = []
    for trace in plate:
        well = trace.name
        trace = pad_trace(trace)
        output_data = dispense_map_chunk.loc[dispense_map_chunk['DestWell'] == well,:].copy().reset_index()
        assert (lo:=len(output_data)) == 0 or lo == 1
        if lo == 0: # empty
            output_data = pd.Series([None] * len(output_data.columns),
                                    index=output_data.columns,
                                    dtype='float64',
                                    )
        if lo == 1: # make df
            output_data = output_data.loc[0,:]

        output_data['Well'] = well
        assert 'Well' in output_data.index
        for i in kwargs:
            output_data[i] = kwargs[i]

        output_data = output_data.append(trace).to_frame().T
        bag_of_df.append(output_data)
    df = pd.concat(bag_of_df)
    return df


def main(config_path):
    output_path = sys.stdout
    assert os.path.exists(config_path)
    if os.path.isdir(config_path):
        files = find(config_path)
        config_path = grep('config.yml', files)
        assert len(config_path) == 1
        config_path = config_path[0]
    config = Config(config_path)
    dispense_map = get_dispense_map(config.echo['picklist'][0], 
                                    config.echo['transfers'])

    HEADERDONE = False
    for plate_n in tqdm(sorted(config.platereader)):
        dispense_map_chunk = dispense_map.loc[dispense_map['dest_plate'] == plate_n, :]
        plate_config = config.platereader[plate_n] #  {'test':{.. , 'control':{..
        test_plate_config = plate_config['test']
        test_plate = PlateData(test_plate_config['path'])
        time, date = test_plate.timestamp
        o = proc_plate(plate=test_plate,
                       output_path=output_path,
                       dispense_map_chunk=dispense_map_chunk,
                       protein=config.protein,
                       append=True,
                       **test_plate_config,
                       )
        if HEADERDONE:
            o.to_csv(output_path, 
                     mode='a', 
                     header=False, 
                     index=False)
        else:
            o.to_csv(output_path, 
                     index=False)
            HEADERDONE = True

        if 'control' in plate_config:
            control_plate_config = plate_config['control']
            control_plate = PlateData(control_plate_config['path'])
            time, date = control_plate.timestamp
            o = proc_plate(plate=control_plate,
                           output_path=output_path,
                           dispense_map_chunk=dispense_map_chunk,
                           protein=None,
                           append=True,
                           **test_plate_config,
                           )
            if HEADERDONE:
                o.to_csv(output_path, 
                         mode='a', 
                         header=False, 
                         index=False)
            else:
                o.to_csv(output_path, 
                         index=False)
                HEADERDONE = True
            #o.to_csv(output_path, 
            #         mode='a', 
            #         header=False, 
            #         index=False)

if __name__ == '__main__':
    main(sys.argv[1])
