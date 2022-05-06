import os
import re
import datetime
from ast import literal_eval
from string import ascii_lowercase, ascii_uppercase
from tqdm import tqdm
#from tqdm.notebook import tqdm

import numpy as np
from scipy.ndimage import gaussian_filter1d, convolve1d
import pandas as pd

import matplotlib.pyplot as plt
plt.style.use('dark_background')
import seaborn as sns

class PlateData:
    def __init__(self,
                 path,
                 **kwargs,
                 ):
        self.__dict__ == {**self.__dict__, **kwargs}
        self.path = path
        self._df = None
    def __len__(self):
        return len(self.df)
    def __iter__(self):
        for i in self.df.index:
            yield self.df.loc[i,:]
    def __repr__(self): # todo - map empty wells
        wells = self.df.index
        p = '\n'.join([i + ' : ' * 24 for i in ascii_uppercase[:16]])
        cols = ''.join(map(lambda i : f' {i} ' if len(str(i)) == 1 else f' {i}',
                       range(1,25)))
        return f'{self.path}\n{self.metadata}\n{cols}\n{p}'
    def __getitem__(self, idx):
        if isinstance(idx, str):
            assert idx in self.df.index
            return self.df.loc[idx, :]
        elif isinstance(idx, list):
            return pd.concat([self[i] for i in idx],
                             axis=1,
                             ).T
        elif isinstance(idx, (slice, int)):
            keys = list(self.df.index)[idx]
            if isinstance(keys, list):
                return pd.concat([self.df.loc[i, :] for i in keys],
                                 axis=1,
                                 ).T
            else:
                return self.df.loc[keys, :]
        else: 
            raise Warning(f'PlateData.__getitem__: issue parsind index {idx}')
    @property
    def df(self):
        if self._df is None:
            self._df = parse(self.path)
        return self._df
    @property
    def metadata(self):
        with open(self.path) as f:
            return ''.join(f.readlines(9))
    @property
    def timestamp(self):
        timestamp_i = lambda s : datetime.datetime.strptime(s.replace('.CSV',''),
                                                            '%d%m%Y,%H%M%S')
        timestamp_j = lambda s : datetime.datetime.strptime(s.split('_')[0].replace('.CSV',''),
                                                            '%d%m%Y,%H%M%S')
        s = os.path.basename(self.path)
        if '_' in s:
            return timestamp_j(s)
        else:
            return timestamp_i(s)

class Cpd:
    ''' Compound : well mapping from picklist
    '''
    def __init__(self,
                 *,
                 df_slice=None,
                 csv=None,
                 search=None,
                 col=None,
                 ):
        if df_slice is not None:
            self.df = df_slice
        else:
            assert csv is not None and search is not None
            assert os.path.exists(csv)
            try:
                df = pd.read_csv(csv)
            except Exception as e:
                print(e)
                raise Warning(f'error reading {csv}')
            matches = []
            for i in df.columns:
                if df[i].dtype == 'O':
                    match_ = df.loc[df[i].str.contains(search),:]
            #matches = [df.loc[i,:] for i in df.index \
            #           if sum(df.loc[i,:].str.contains(search)) > 0 ]
            if len(matches) > 0:
                self.df = pd.concat(matches)
            else:
                raise Warning(f'No results in {csv} for {search}')
    def __len__(self):
        return len(self.df)
    def __iter__(self):
        for i in self.df.index:
            yield self.df.loc[i, :]
    def __repr__(self): 
        pass
    def __getitem__(self, idx):
        assert idx in self.df.index, f'{idx} not in index: {self.df.index}'
        return self.df.loc[idx,:]
    @property
    def metadata(self):
        pass
    

class Screen:
    ''' class for handling a set of plates 
        and their picklist, mapping wells to
        compounds and volumes.
    '''
    def __init__(self,
                 *,
                 picklist,
                 files,
                 exceptions=None,
                ):
        self.files = files
        def _read_picklist(picklist):
            if isinstance(picklist, str):
                if os.path.exists(picklist):
                    picklist = pd.read_csv(picklist, index_col=0)
                    return picklist
                else:
                    raise Exception(f'OS: picklist {picklist} does not exist ')
            elif isinstance(picklist, pd.DataFrame):
                return picklist

        if isinstance(picklist, (list, tuple)):
            raise Warning('Read multiple picklist not implemented')
            #self.picklist = {i:_read_picklist(i) for i in picklist}
        elif isinstance(picklist, str):
            # assumes that one picklist applies to both
            self.picklist = _read_picklist(picklist) 
            
        cpd_num = lambda s : re.search('(S[0-9]+)', s).groups()[0]
        # dict of {cpd:pd.DataFrame, ... } - picklist wells
        self.data  = {i:PlateData(j) for i, j in enumerate(files)}
        self.cpds = {cpd_num(i):{'picklist':self.picklist.loc[self.picklist['Cpd'] == i, :],
                                 'wells':None}
                             for i in self.picklist['Cpd'].unique()}
        
    def __len__(self):
        return len(self.cpds)
    def __getitem__(self, idx):
        if isinstance(idx, str):
            assert idx in self.cpds.keys()
            return self.proc(idx)
        elif isinstance(idx, (slice, float)):
            keys = list(self.data.keys())[idx]
            return [self.data[i] for i in keys]
    def __iter__(self):
        pass
    def __repr__(self):
        return f'Screen\n{self.__len__()} compounds\n{len(self.files)} plate reader data files'
    @property
    def cpdNames(self):
        return list(self.cpds.keys())
    def proc(self, cpd):
        data_ = self.cpds[cpd]
        picklist = data_['picklist']
        plate_name_ = list(set(picklist['Destination Plate Name']))
        assert len(plate_name_) == 1, 'multiple files look up not implemented'
        plate_name = plate_name_[0]
        wells = picklist['DestWell'].to_list()
        vols = picklist['Transfer Volume /nl'].to_list()
        data = None
        #return {'plateName': plate_name,
        #        'wells'    : wells,
        #        'vols'     : vols}


def parse(path):
    def _proc(df):
        if 'Unnamed: 0' in df.columns or 'Unnamed: 1' in df.columns:
            wells = [f'{i}{j}' for i,j in zip(df['Unnamed: 0'], 
                                              df['Unnamed: 1'])]
            df.index =  wells
            for i in ['Unnamed: 0',
                      'Unnamed: 1',
                      'Wavelength [nm]',
                      'Content']:
                if i in df.columns:
                    df.drop(i,axis=1, inplace=True)
        df.dropna(axis=1, inplace=True)
        df=df.replace('overflow', 3.5)
        df.columns = list(map(int, df.columns)) # wavelengths
        return df
    try:
        df = pd.read_csv(path, skiprows=6)
        return _proc(df)
    except:
        df = pd.read_csv(path, skiprows=7)
        return _proc(df)

def timestamp(s):
    timestamp_i = lambda s : datetime.datetime.strptime(s.replace('.CSV',''),
                                                        '%d%m%Y,%H%M%S')
    timestamp_j = lambda s : datetime.datetime.strptime(s.split('_')[0].replace('.CSV',''),
                                                        '%d%m%Y,%H%M%S')
    if '_' in s:
        return timestamp_j(s)
    else:
        return timestamp_i(s)

def proc(df):
    cols = df.columns
    idx = df.index
    dfZ800 = df.sub(df.loc[:,800], axis=0)
    dfSmth = pd.DataFrame(gaussian_filter1d(dfZ800, 3),
                          index=idx,
                          columns=cols)
    return dfSmth


def plotTraces(x,     # df
               save_path=None,
               cpd=None,   # name
               vols=None,  #
               save=False,
               size=(12,8),
               **kwargs,
               ):
    if save_path is not None:
        if not os.path.exists('img'):
            os.mkdir('img')
    plt.figure(figsize=size)
    if vols is not None:
        for row_, vol_ in zip(x.index, vols):
            plt.plot(x.loc[row_,:], 
                     c=plt.cm.cool(vol_/2000), 
                     label=f'{row_} {vol_} uM')
    else:
        for row_ in x.index:
            plt.plot(x.loc[row_,:],
                     alpha=0.5,
                     lw=1,
                     label=f'{row_}',
                     **kwargs,
                     )
    plt.xlim(280,800)
    plt.ylim(0, 1)
    if cpd is not None:
        plt.title(cpd)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Absorbance')
    if len(x.index) < 10:
        plt.legend()
    if save or save_path is not None:
        if save_path is not None:
            plt.savefig(os.path.join('img', f'{save_path}.png'))
        else:
            plt.savefig(os.path.join('img', f'{cpd}.png'))
        plt.close()
    else:
        plt.show()
    

def plot_set(wells, title=''):
    from tqdm import tqdm
    assert len(wells) > 0, 'empty set given'
    if len(wells) < 8:
        fig, ax = plt.subplots(len(wells), 
                               1, 
                               figsize=(6,16))
    else:
        nRows=8
        fig, ax = plt.subplots(len(wells)//nRows, 
                               nRows, 
                               figsize=(16,128*(len(wells) / len(wells))))
        #fig, ax = plt.subplots(len(wells)//(nRows:=8), 
        #                       nRows, 
        #                       figsize=(16,128*(len(wells) / len(wells))))

    for cpd, ax_ in zip(wells, ax.flatten()):
        _plate = list(set(wells[cpd]['Destination Plate Name']))[0]
        _wells = wells[cpd]['DestWell'].to_list()
        _vols = wells[cpd]['Transfer Volume /nl'].to_list()
        data_ = data[next(filter(lambda item : data[item]['pname'] == _plate, data))]
        data_cpd = data_['x'].loc[_wells,:]

        x = proc(data_cpd) 

        for row_, vol_ in zip(x.index, _vols):
            ax_.plot(x.loc[row_,:], c=plt.cm.cool(vol_/2000))

        ax_.set_xlim(280,800)
        ax_.set_title(cpd)
        ax_.set_xlabel('Wavelength (nm)')
        ax_.axis('off')

    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_report(*,
                traces=None,
                diff=None,
                name=None,
                smiles=None,
                x=None,
                y=None,
                mm=None,
                save_path=None,
                **kwargs,
                ):
    fig, ax = plt.subplots(2,2, figsize=(6,6))
    if traces is not None:
        if x is None:
            x = dict(zip(traces.index, range(len(traces)))) # problem later it .T
        p = ax[0,0]
        for trace_, conc_, well_ in zip(traces.index, 
                                        x.values(),
                                        x.keys()):
            p.plot(traces.loc[well_,:], 
                   c=plt.cm.cool(conc_/max(x.values())))
        p.set_xlim(280,800)
        p.set_ylim(-0.1,0.25)
        p.set_title(name)
        p.set_xlabel('Wavelength (nm)')
        #p.axis('off')
    if diff is not None:
        p = ax[0,1]
        for diff_, conc_, well_ in zip(diff.index, 
                                       x.values(),
                                       x.keys()):
            p.plot(diff.loc[well_,:], 
                   c=plt.cm.cool(conc_/max(x.values())),
                   label=f"{conc_} uM")
        p.legend()
        p.set_xlim(280,max(diff.columns))
        p.set_ylim(-0.1,0.2)
        p.set_title(name)
        p.set_xlabel('Wavelength (nm)')
        #p.axis('off')
    if x is not None and y is not None:
        p=ax[1,0]
        p.scatter(x.values(), y.values())
        if mm is not None:
            x_ = np.linspace(min(x.values()), 
                             max(x.values()), 
                             32)
            from pipeline import michaelis_menten
            y_ = michaelis_menten(x_, mm['km'], mm['vmax'])
    if smiles is not None:
        assert isinstance(smiles, str)
        from rdkit import Chem
        from rdkit.Chem import Draw
        p=ax[1,1]
        m = Chem.MolFromSmiles(smiles)
        im = Draw.MolToImage(m)
        p.imshow(im)
        p.axis('off')

    plt.tight_layout()

    if save_path is not None:
        directory = os.path.dirname(save_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        plt.savefig(save_path)

    plt.close()
