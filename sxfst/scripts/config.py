#!/usr/bin/env python
import sys
import os
import re
import yaml
from pprint import pprint
import argparse

'''
generates a config.yml for an experiment, limited reliability because of
my incosnsistent filing methods.
'''

NPLATES_PER_RUN = 15

def contains(regex, string, **kwargs):
    return re.search(regex, string, **kwargs) is not None

def grep(regex, data, *args):
    # string or list
    if isinstance(data, str):
        search = re.search(regex, data, *args)
        if search is not None:
            if len(search.groups()) == 1:
                return search.group(1) # extract
            else:
                return search.group()
    elif isinstance(data, list):
        return [i for i in data if contains(regex, i)]
    else:

        raise Warning()

def find(path):
    assert os.path.exists(path)
    return [i for i in \
            os.popen(f'find {path}').read().replace('//','/').split('\n') \
            if i != '']

def filter_files(paths):
    return [i for i in paths if os.path.isfile(i)]

def get_plate_metadata(path):
    assert os.path.exists(path)
    assert os.path.isfile(path)
    with open(path) as f:
        head = f.read(1024).replace('\\','/')
    test_run_nos = grep('Test run no.: (\d{4})', head)
    time = grep('Time:\s(.*),', head)
    time = grep('[0-9]{2}:[0-9]{2}:[0-9]{2}', head)
    date = grep('Date:\s*(.*),', head)
    user = grep('User:\s*([A-Za-z]*),', head)
    id1 = grep('ID1:\s*([A-Za-z0-9]*),', head)
    id2 = grep('ID2:\s*([A-Za-z0-9]*),', head)
    machine = grep('(BMG/[A-Za-z0-9]*)/', head).replace('/', ' ')

    return {'test_run_no':test_run_nos,
            'date':date,
            'time':time,
            'id1':id1,
            'id2':id2,
            'machine':machine,
            'user':user,
            'path':path}


def main(args):
    files = find(args.input)
    cfg = {}

    echo_files = [i for i in grep('echo', files) if os.path.isfile(i)]
    cfg['echo'] = {'picklist':grep('picklist', echo_files),
                   'protocol':grep('.*epr$', echo_files),
                   'surveys':grep('Survey', echo_files),
                   'transfers':grep('Transfer', echo_files)}
    uv_files = [i for i in grep('uv-vis', files) if os.path.isfile(i)]
    cfg['uv-vis'] = {'pre-dilution spec':grep('pre',uv_files),
                     'post-dilution spec':grep('post',uv_files),
                     }
    cfg['nb'] = [i for i in grep('ipynb', files) if os.path.isfile(i)]


    platereader_files = [i for i in grep('platereader', files) \
                        if os.path.isfile(i)]
    plates_ = {(a:=get_plate_metadata(i))['test_run_no']:a \
            for i in platereader_files}
    plates = {i:plates_[i] for i in sorted(plates_)} # run order
    if args.controls is not None:
        assert os.path.exists(args.controls)
        if os.path.isdir(args.controls):
            ctrl_files = find(args.controls)
            ctrl_config_path = grep('config.yml', ctrl_files)
        elif os.path.isfile(args.controls):
            ctrl_config_path = args.controls
        with open(ctrl_config_path) as f:
            ctrl_config = yaml.full_load(f)
        if 'platereader' not in cfg.keys():
            cfg['platereader'] = {}
        for i, j in enumerate(plates, 1):
            plate_name = f'plate_{i}'
            test_plate_info = plates[j]
            ctrl_plate_info = ctrl_config['platereader'][plate_name]['control']
            ctrl_plate_rel_path = os.path.join(os.path.dirname(args.controls), 
                                               ctrl_plate_info['path'])
            ctrl_plate_info['path'] = ctrl_plate_rel_path
            cfg['platereader'][plate_name] = {'test':plates[j],
                                                'control':ctrl_plate_info}
    elif len(plates) == NPLATES_PER_RUN:
        # no controls
        # assumes read order was plate order
        # should be ok for my experiments
        cfg['platereader'] = {f'plate_{i}':{'test':plates[j],
                        'control':None,} for i,j in enumerate(plates, 1)}
    elif len(plates) > NPLATES_PER_RUN:
        if len(plates) % NPLATES_PER_RUN == 0:
            # control or extra set
            if len(plates) / NPLATES_PER_RUN == 2:
                # assume control first
                ctrl_keys = (a:=sorted(plates.keys()))[:NPLATES_PER_RUN]
                test_keys = a[NPLATES_PER_RUN:]
                cfg['platereader'] = {f'plate_{i}':{'test':plates[k],
                                                    'control':plates[j],
                                                    }\
                    for i,(j,k) in enumerate(zip(ctrl_keys, test_keys), 1)}
            elif len(plates) / NPLATES_PER_RUN == 3:
                # assume control -> unclassified protein -> good protein
                ctrl_keys = (a:=sorted(plates.keys()))[:NPLATES_PER_RUN]
                unclassified_keys = a[NPLATES_PER_RUN:NPLATES_PER_RUN*2]
                test_keys = a[NPLATES_PER_RUN*2:]
                cfg['platereader'] = {f'plate_{i}':{'test':plates[k],
                                                    'control':plates[j],
                                                    'unclassified':plates[l],
                                                    }\
                    for i,(j,k,l) in enumerate(zip(ctrl_keys, 
                                                   test_keys, 
                                                   unclassified_keys), 
                                    1)}
            else:
                raise Warning(f'Not impelented for {len(plates)}, sorry')

    print(yaml.dump(cfg))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', 
                        help='directory path of experiment')
    parser.add_argument('-c', '--controls',
                        help='the config that contains the controls')
    args = parser.parse_args()
    main(args)
    #main(sys.argv[1:])
