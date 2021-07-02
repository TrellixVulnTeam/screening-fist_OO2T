import os
import os.path as osp
import ast
from pprint import pprint
import pandas as pd
from tqdm import tqdm

import plates
import echo


def get_planned_vols(map_df):
    return {i:ast.literal_eval(j) for i,j in zip(map_df['well'], map_df['contents'])}

def get_actual_vols(exceptions, map_df, plate_name):
    exceptions_chunk = exceptions.loc[exceptions['Destination Plate Name'] == plate_name,:]
    planned_vols = get_planned_vols(map_df)
    for i,j in zip(exceptions_chunk['Destination Well'], exceptions_chunk['Actual Volume']):
        well = planned_vols[i]
        key = [i for i in well.keys() if i != 'dmso'][0]
        well[key] = j
    return planned_vols

def actual_vols_to_df(actual_vols):
    values = actual_vols.values()
    def get_contents(d):
        if 'dmso' in d.keys():
            dmso = d['dmso']
            if len(d.keys()) > 1:
                cpd_name = [i for i in d.keys() if i !='dmso'][0]
                cpd_vol = d[cpd_name]
            else:
                cpd_name = None
                cpd_vol = 0
        else:
            dmso = 0
            if len(d.keys()) == 1: # all compound
                cpd_name = next(iter(d.keys()))
                cpd_vol =  next(iter(d.values()))
            else: # empty wells 
                cpd_name = None
                cpd_vol = 0 
        return dmso, cpd_vol, cpd_name
    return pd.DataFrame([get_contents(i) for i in actual_vols.values()],
                       index = actual_vols.keys(),
                       columns = ['dmso vol', 'cpd vol', 'cpd name'])


def process(df,exceptions, name):
    x = actual_vols_to_df(get_actual_vols(exceptions,df,name))
    total_vol = x['dmso vol'] + x['cpd vol'] + 28
    x['cpd conc µM'] = x['cpd vol'] * 10_000 / total_vol# (v1*c1) / v2
    return x


def main():
    # picklists
    dmso_picklist = pd.read_csv('../design/dmso-src-picklist.csv', index_col=0)
    cpd_picklist = pd.read_csv('../design/cpd-src-picklist.csv', index_col=0)
    # plate maps
    bsa_blank_map = pd.read_csv('../design/bsa-blank-map.csv', index_col=0)
    bsa_test_map = pd.read_csv('../design/bsa-test-map.csv', index_col=0)
    kpi_blank_map = pd.read_csv('../design/ctrl-blank-map.csv', index_col=0)
    kpi_test_map = pd.read_csv('../design/ctrl-test-map.csv', index_col=0)
    exceptions = pd.read_csv('../data/20210608_src-except.csv')

    name_map = {}
    for i,j in zip(exceptions['Destination Plate Name'].unique(),
                  cpd_picklist['Destination Plate Name'].unique()):
        print(i,'----',j)
        name_map[i]=j


    bsa_blank_map_actual = process(bsa_blank_map, exceptions,'Destination[2]')
    bsa_test_map_actual = process(bsa_test_map, exceptions,'Destination[3]')
    kpi_blank_map_actual = process(kpi_blank_map, exceptions,'Destination[4]')
    kpi_test_map_actual = process(kpi_test_map, exceptions,'Destination[5]')


    bsa_blank_map_actual.to_csv('bsa_blank_map_actual.csv')
    bsa_test_map_actual.to_csv('bsa_test_map_actual.csv')
    kpi_blank_map_actual.to_csv('kpi_blank_map_actual.csv')
    kpi_test_map_actual.to_csv('kpi_test_map_actual.csv')


    plates_data = {'bsa-blank':{0:'../data/plates/kpi-bsa-blank-t0.CSV',
                                1:'../data/plates/bsa-blank-t1.CSV',
                                2:'../data/plates/bsa-blank-t2.CSV',
                                3:'../data/plates/bsa-blank-t3.CSV'},
                  'bsa-bm3':{0:'../data/plates/kpi-bsa-bm3-t0.CSV',
                             1:'../data/plates/bsa-bm3-t1.CSV',
                             2:'../data/plates/bsa-bm3-t2.CSV',
                             3:'../data/plates/bsa-bm3-t3.CSV'},
                   'kpi-blank':{0:'../data/plates/kpi-blank-t0.CSV',
                             1:'../data/plates/kpi-blank-t1.CSV',
                             2:'../data/plates/kpi-blank-t2.CSV',
                             3:'../data/plates/kpi-blank-t3.CSV'},
                   'kpi-bm3':{0:'../data/plates/kpi-bm3-t0.CSV',
                             1:'../data/plates/kpi-bm3-t1.CSV',
                             2:'../data/plates/kpi-bm3-t2.CSV',
                             3:'../data/plates/kpi-bm3-t3.CSV'}
                  }

    #test paths
    plate_objs = {}
    for i in plates_data:
        plate_objs[i] = {}
        for j in plates_data[i]:
            plate_objs[i][j] = plates.uv384.UV384m4(plates_data[i][j]) 

    keys = {'kpi':['kpi-bm3','kpi-blank'],'bsa':['bsa-bm3','bsa-blank']}


    if not osp.exists('reports-notnorm'):
        os.makedirs('reports-notnorm', exist_ok = True)

    for test, ctrl in keys.values():
        for timepoint, (test_plate, ctrl_plate) in enumerate(zip(plate_objs[test].values(), plate_objs[ctrl].values())):
            for block_num, (test_block, ctrl_block) in tqdm(enumerate(zip(test_plate.blocks, ctrl_plate.blocks)), total = len(test_plate.blocks)):
                test_concs = kpi_test_map_actual.loc[test_block.index,'cpd conc µM']
                ctrl_concs = kpi_blank_map_actual.loc[ctrl_block.index,'cpd conc µM']
                block = plates.Block(data=test_block, 
                                     control=ctrl_block, 
                                     test_concs=test_concs,
                                     ctrl_concs = ctrl_concs)
                plates.report(block, save_path = osp.join('reports-notnorm',f'{test}-block{block_num}-t{timepoint}.png'))

    if not osp.exists('reports-smooth-norm'):
        os.makedirs('reports-smooth-norm', exist_ok = True)

    for test, ctrl in keys.values():
        for timepoint, (test_plate, ctrl_plate) in enumerate(zip(plate_objs[test].values(), plate_objs[ctrl].values())):
            for block_num, (test_block, ctrl_block) in tqdm(enumerate(zip(test_plate.blocks, ctrl_plate.blocks)), total = len(test_plate.blocks)):
                test_concs = kpi_test_map_actual.loc[test_block.index,'cpd conc µM']
                ctrl_concs = kpi_blank_map_actual.loc[ctrl_block.index,'cpd conc µM']
                block = plates.Block(data=test_block, 
                                     control=ctrl_block, 
                                     test_concs=test_concs,
                                     ctrl_concs = ctrl_concs,
                                     norm410 = True,
                                     smooth = True)
                plates.report(block, save_path = osp.join('reports-smooth-norm',f'{test}-block{block_num}-t{timepoint}.png'))

    # todo
    # why are some traces empty?
    # is interpolation working properly?
    # test 2 factor scaling
    # mm fit: bootstrap with replacement

if __name__ == '__main__':
    main()
