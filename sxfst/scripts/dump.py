import heapq
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters
from tqdm import tqdm

def filter_HL(smiles):
    # filter smiles by herbicide likeness rules Hao
    # todo: n aromatic bonds < 17
    mw_cutoff = 435
    logp_cutoff = 6
    hba_cutoff = 6
    hbd_cutoff = 2
    rotatable_cutoff =  9
    n_aromatic_bonds_cutoff = 17
    
    makemols = lambda smiles : Chem.AddHs(Chem.MolFromSmiles(smiles))
    n_aromatic_bonds = lambda m : sum([b.GetIsAromatic() for b in m.GetBonds()])
    mols = [makemols(i) for i in smiles]
    props = {s:{'mw':Chem.CalcExactMolWt(i), 
                'logp': Crippen.MolLogP(i),
                'hba': Chem.CalcNumLipinskiHBA(i),
                'hbd': Chem.CalcNumLipinskiHBD(i),
                'rot': Chem.CalcNumRotatableBonds(i),
                'aroB': n_aromatic_bonds(i)}
                for i,s in zip(mols, smiles)}
    
    prop_filter = lambda s : props[s]['mw'] <= mw_cutoff \
                        and props[s]['logp'] <= logp_cutoff \
                        and props[s]['hba'] <= hba_cutoff\
                        and props[s]['hbd'] <= hbd_cutoff\
                        and props[s]['rot'] <= rotatable_cutoff\
                        and props[s]['aroB'] <= n_aromatic_bonds_cutoff

    return [i for i in smiles if prop_filter(i)]


def pick(smiles, n):
    picker = SimDivFilters.MaxMinPicker()
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in mols]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    idx = picker.LazyPick(fn, len(smiles), n)
    return [smiles[i] for i in idx]


def pick_set(filteredSMILES, plates, n):
    selectionSMILES = pick(filteredSMILES, n)
    selection = pd.concat([plates.loc[plates['CanonicalSMILES'] == i, :] for i in selectionSMILES])
    locations = selection[['Plate','Well','Product Name']]
    return locations.sort_values(['Plate','Well'])

def main():
    smiles = pd.read_csv('lib/fda-canonicalSMILES-deduplicated.csv', index_col=-1)['CanonicalSMILES']
    

    plates = pd.read_csv('lib/fda.csv')
    in_stock = ['HY-L022-1', 'HY-L022-2']
    lib = plates.loc[plates['Plate'].isin(in_stock), :]

    canonicalSMILES = []
    for i in lib['Product Name']:
        if i not in smiles.index:
            canonicalSMILES.append(None)
        else:
            canonicalSMILES.append(smiles.loc[i])

    lib = pd.concat([lib, pd.Series(canonicalSMILES, name = 'CanonicalSMILES')], axis = 1, join = 'inner')
    lib.drop(([i for i, j in zip(lib.index, lib['CanonicalSMILES']) if j is None]), inplace = True)


    filtered = filter_HL(lib['CanonicalSMILES'])
    selection = pick_set(filtered, lib, 48)
    print(selection)
    selection.to_csv('selection.csv', index=False)

if __name__ == '__main__':
    main()
import pandas as pd

def main():
    df = pd.read_excel('selleckchem-plate.xlsx',
            sheet_name = 'L1300-FDA-978cpds')
    df = df.loc[:,['Item Name','CatalogNumber', 'SMILES', 'Rack Number', 'Plate Location']]
    df.to_csv('layouts.csv')

if __name__ == '__main__':
    main()
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters


def filter_HL(smiles):
    # filter smiles by herbicide likeness rules Hao
    # todo: n aromatic bonds < 17
    mw_cutoff_min = 100
    mw_cutoff_max = 435
    logp_cutoff = 6
    hba_cutoff = 6
    hbd_cutoff = 2
    rotatable_cutoff =  9
    n_aromatic_bonds_cutoff = 17
    
    makemols = lambda smiles : Chem.AddHs(Chem.MolFromSmiles(smiles))
    n_aromatic_bonds = lambda m : sum([b.GetIsAromatic() for b in m.GetBonds()])
    mols = [makemols(i) for i in smiles]
    props = {s:{'mw':Chem.CalcExactMolWt(i), 
                'logp': Crippen.MolLogP(i),
                'hba': Chem.CalcNumLipinskiHBA(i),
                'hbd': Chem.CalcNumLipinskiHBD(i),
                'rot': Chem.CalcNumRotatableBonds(i),
                'aroB': n_aromatic_bonds(i)}
                for i,s in zip(mols, smiles)}
    
    prop_filter = lambda s : props[s]['mw'] <= mw_cutoff_max \
                        and props[s]['mw'] >= mw_cutoff_min \
                        and props[s]['logp'] <= logp_cutoff \
                        and props[s]['hba'] <= hba_cutoff\
                        and props[s]['hbd'] <= hbd_cutoff\
                        and props[s]['rot'] <= rotatable_cutoff\
                        and props[s]['aroB'] <= n_aromatic_bonds_cutoff
    return [i for i in smiles if prop_filter(i)]

def pick(smiles, n):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in mols]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    picker = SimDivFilters.MaxMinPicker()
    idx = picker.LazyPick(fn, len(smiles), n)
    return [smiles[i] for i in idx]

def lookup(smiles, df):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main():
    df = pd.read_csv('layouts.csv')
    herbicideLike = filter_HL(df['SMILES'])
    selection = lookup(pick(herbicideLike, 39), df)
    selection.drop('Unnamed: 0', axis = 1, inplace=True) 
    selection.to_csv('herbicide-like-selection.csv', index=False)

if __name__ == '__main__':
    main()
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters


def filter_HL(smiles):
    # filter smiles by herbicide likeness rules Hao
    mw_cutoff_min = 100
    mw_cutoff_max = 435
    logp_cutoff = 6
    hba_cutoff = 6
    hbd_cutoff = 2
    rotatable_cutoff =  9
    n_aromatic_bonds_cutoff = 17
    
    makemols = lambda smiles : Chem.AddHs(Chem.MolFromSmiles(smiles))
    n_aromatic_bonds = lambda m : sum([b.GetIsAromatic() for b in m.GetBonds()])
    mols = [makemols(i) for i in smiles]
    props = {s:{'mw':Chem.CalcExactMolWt(i), 
                'logp': Crippen.MolLogP(i),
                'hba': Chem.CalcNumLipinskiHBA(i),
                'hbd': Chem.CalcNumLipinskiHBD(i),
                'rot': Chem.CalcNumRotatableBonds(i),
                'aroB': n_aromatic_bonds(i)}
                for i,s in zip(mols, smiles)}
    
    prop_filter = lambda s : props[s]['mw'] <= mw_cutoff_max \
                        and props[s]['mw'] >= mw_cutoff_min \
                        and props[s]['logp'] <= logp_cutoff \
                        and props[s]['hba'] <= hba_cutoff\
                        and props[s]['hbd'] <= hbd_cutoff\
                        and props[s]['rot'] <= rotatable_cutoff\
                        and props[s]['aroB'] <= n_aromatic_bonds_cutoff
    return [i for i in smiles if prop_filter(i)]

def lookup(smiles, df):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main():
    df = pd.read_csv('layouts.csv')
    herbicideLike = filter_HL(df['SMILES'])
    #herbicideLike = lookup(herbicideLike, df)
    print(len(herbicideLike))
    # herbicideLike.to_csv('herbicide-like.csv')

if __name__ == '__main__':
    main()
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit import SimDivFilters

def pick(smiles, n):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in smiles]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    picker = SimDivFilters.MaxMinPicker()
    idx = picker.LazyPick(fn, len(smiles), n)
    return pd.Series([smiles[i] for i in idx], index = idx)


def main():
    df = pd.read_csv('fda.csv') # smiles errors, scrape cas numbers
    print(pick(df['SMILES'], 10))

if __name__ == '__main__':
    main()
import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def process(df):
    df.index = df.iloc[:,1]
    df = df.iloc[:,2:]
    return df

def plot(df, name = None):
    plt.figure(figsize=(15,5))
    sns.heatmap(np.random.randn(8,12), 
            annot=df, 
            fmt = '',
            cbar = False,
            annot_kws={'fontsize':10,
                'rotation':30})
    plt.xticks(range(1,13), range(1,13))
    plt.yticks(range(8), list('abcdefgh'), rotation = 0)
    plt.title(name)
    if name != None:
        plt.savefig(name + '.png')
    else:
        plt.show()



def main():
    csvs = [i for i in os.listdir() if 'hy-' in i and 'csv' in i]
    dfs = [pd.read_csv(i, header=None) for i in csvs]
    dfs = [process(i) for i in dfs]
    for i, j in tqdm(zip(dfs, csvs)):
        try:
            plot(i, name = j.split('.')[0])
        except  Exception as e:
            print(e)

if __name__ == '__main__':
    main()
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

def main():
    df = pd.read_csv('fda-canonical-smiles.csv')[['name','CanonicalSMILES']].drop_duplicates()

    repeats = df['name'].value_counts().loc[df['name'].value_counts() > 1]

    mols = []
    names = []
    for i in repeats.index:
        x = df.loc[df['name'] == i, :]
        for j in x['CanonicalSMILES']:
            mols.append(Chem.MolFromSmiles(j))
            names.append(i)
    im = Draw.MolsToGridImage(mols, legends = names)
    im.save('duplicated-smiles.png')

if __name__ == '__main__':
    main()
import time
from io import StringIO
import requests 
import pandas as pd
from tqdm import tqdm

def get_cid(cas):
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=pccompound&UID={cas}&rettype=SMILES'
    r = requests.get(url)
    print(r.status_code)
    if r.status_code == 200:
        print(r.text)
    if r.status_code == 414:
        print('url too long')

def search_name(name):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/CSV'
    r = requests.get(url)
    if r.status_code == 200:
        data = r.text
        df = pd.read_csv(StringIO(data))
        df['name'] = name
        return df
    else:
        print(name, r.status_code)
        pd.DataFrame([], columns = ['CID','CanonicalSMILES','name'])


def get_smiles(cid):
    pass

def main():
    df = pd.read_csv('fda.csv')
    cas = df['CAS Number']
    names = df['Product Name']


    results = pd.DataFrame([], columns = ['CID','CanonicalSMILES','name'])
    for i in tqdm(names):
        try:
            results = results.append(search_name(i))
        except :
            time.sleep(60)
            results = results.append(search_name(i))
    results.to_csv('fda-canonical-smiles.csv')

if __name__ == '__main__':
    main()
import string
import numpy as np
import pandas as pd

class SourcePlateCompound:
    def __init__(self, coumpound_name, wells, well_vol = None, ldv = True):
        # give vol in µl
        assert isinstance(wells, list)
        self.compound = coumpound_name
        self.ldv = ldv
        if self.ldv:
            self.MaxWellVol = 12 * 1000 #nl
            self.MinWellVol = 2.5 * 1000 #nl + safety
        else:
            self.MaxWellVol = 65 * 1000 #nl
            self.MinWellVol = 15 * 1000 #nl + safety

        if well_vol is None:
            self.well_vol = self.MaxWellVol
        else:
            self.well_vol = well_vol * 1000 # µl -> nl
        assert self.well_vol <= self.MaxWellVol
        self.wells = self.FillWells(wells)

    def FillWells(self,wells):
        output = {}
        for i in wells:
            output[i] = self.well_vol
        return output

    @property
    def AvailableVolume(self):
        return sum(self.wells.values())

    def Sample(self,vol):
        if vol %2.5 !=0:
            print('Transfer vol not a multiple of 2.5')
        sample = {}
        # take what you can then move on to the next well
        for well in self.wells:
            if self.wells[well] > self.MinWellVol:
                if vol < self.wells[well] - self.MinWellVol:
                    self.wells[well] -= vol
                    sample[well] = vol
                    vol -=vol
                    break
                else:
                    AvailableVol = self.wells[well]-self.MinWellVol
                    sample[well] = AvailableVol
                    self.wells[well] -= AvailableVol
                    vol -= AvailableVol
            else:
                pass # next well
        if vol !=0:
            print(f'{self.compound}: \tVol not reached \t {self.AvailableVolume / 1000} µl')
        return sample

class Block():
    def __init__(self, Compound, DMSO, WorkingVol):
        self.WorkingVol = WorkingVol *1000 # convert to nl
        self.ProteinConc = 10 #ish
        self.K = 3 # prevents duplicates of zero values at 20 ul working vol
        self.Percent_DMSO = 0.05 # as a fraction of 1
        self.Compound = Compound
        self.DMSO = DMSO
        self.TestWells = ['X'+str(i) for i in range(1,9)]
        self.Transfers = self.MakeTransfer()
        
    def MakeTransfer(self):
        compound_vol = np.linspace(0,1,8)**self.K
        compound_vol *= self.Percent_DMSO
        compound_vol *= self.WorkingVol
        compound_vol = 2.5* np.round(compound_vol/2.5)
        TotalDMSOVol = np.round((self.Percent_DMSO * self.WorkingVol)/2.5) *2.5
        DMSO = TotalDMSOVol - compound_vol

        vols = {self.Compound:[i for i in compound_vol],\
             self.DMSO:[i for i in DMSO]}
        
        output = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        for vol_cpd, vol_dmso ,testwell in zip(vols[self.Compound], vols[self.DMSO],self.TestWells):
            cpd_transfer = self.Compound.Sample(vol_cpd)
            dmso_transfer = self.DMSO.Sample(vol_dmso)
            for i,j in zip(cpd_transfer, dmso_transfer):
                temp = pd.DataFrame(\
                [[i,testwell,cpd_transfer[i]],\
                [j,testwell,dmso_transfer[j]]],\
                columns = ['SrcWell','Destination Well','Volume'])
                output = output.append(temp)
        return output.reset_index(drop=True)
            
class AssayPlate():
    def __init__(self):
        self.blocks = {}
        self.TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        self.wells = [f'{i}{j}' for i in string.ascii_uppercase[:16] for j in range(1, 25)]
    def AddBlocks(self,block):
        count = len(self.blocks)+1
        self.blocks[count] = block

    def MapWells(self):
        TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        alphabet = string.ascii_uppercase
        for BlockNumber, wells in zip(self.blocks, 
                [self.wells[i*8:i*8 + 8] for i in range(round(len(self.wells)/8))]):
            mappings = dict(zip(['X'+str(i) for i in range(1,9)], wells))
            transfers = self.blocks[BlockNumber].Transfers
            transfers['Destination Well'] = transfers['Destination Well'].replace(mappings)
            self.TransferPlan = self.TransferPlan.append(transfers)
            self.TransferPlan.reset_index(inplace=True,drop=True)
        self.TransferPlan = self.TransferPlan.loc[self.TransferPlan['Volume'] > 0]
    
def main():
    df = pd.read_csv('herbicide-like-selection.csv')
    dmso_sourcewells = [f'{i}{j}' for i in string.ascii_uppercase[:6] for j in range(1, 25)]

    DMSO = SourcePlateCompound('DMSO',dmso_sourcewells, well_vol = 11, ldv=True)
    cpd_src_wells = [f'{i}{j}' for i in string.ascii_uppercase[:16] for j in range(1, 25)]
    cpd_src_wells = [i for i in cpd_src_wells if i not in dmso_sourcewells]
    split = lambda l, n : [l[i*n:i*n + n] for i in range(round(len(l) / n))]
    cpd_src_tuple = split(cpd_src_wells, 2)
    CPDS = [SourcePlateCompound(i, j, well_vol = 11, ldv = True) for i, j in zip(df['Item Name'], cpd_src_tuple)]

    for repeat in range(4):
        assayplate = AssayPlate()
        for i in CPDS:
            block = Block(i,DMSO,30)
            assayplate.AddBlocks(block)
        assayplate.MapWells()
        transferMap = assayplate.TransferPlan
        transferMap.to_csv(f'picklist{repeat}.csv', index=False)

    '''
    print('DMSO ',DMSO.AvailableVolume / 1000, ' µl')
    for i in CPDS:
        print(i.compound, i.AvailableVolume / 1000, ' µl')
    '''
    src_layout = pd.DataFrame([['DMSO', i, 11] for i in DMSO.wells] +\
            [[i.compound, j, 11] for i in CPDS for j in i.wells], 
            columns = ['well', 'contents', 'vol'])

    print(src_layout)
    for i in CPDS:
        print(i.AvailableVolume/1000, i.wells)
    src_layout.to_csv('src_layout.csv', index=False)


if __name__ == '__main__':
    main()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import string



class SourcePlateCompound():
    def __init__(self,coumpound_name,wells,ldv = True):
        self.compound = coumpound_name
        self.ldv = ldv
        if self.ldv:
            self.MaxWellVol = 12 * 1000 #nl
            self.MinWellVol = 2.5 * 1000 #nl + safety
        else:
            self.MaxWellVol = 65 * 1000 #nl
            self.MinWellVol = 20 * 1000 #nl + safety
        self.wells = self.FillWells(wells)

    def FillWells(self,wells,):
        output = {}
        
        if self.ldv:
            for i in wells:
                output[i] = 11 * 1000 #self.MaxWellVol #ul
        else:
            for i in wells:
                output[i] = 60 * 1000
        return output

    def AvailableVolume(self):
        return sum(self.wells.values())

    def Sample(self,vol):
        if vol %2.5 !=0:
            print('Transfer vol not a multiple of 2.5')
        sample = {}
        # take what you can then move on to the next well
        for well in self.wells:
            if self.wells[well] > (self.MinWellVol + 5000): # extra safety margin
                
                if vol < self.wells[well]:
                    self.wells[well] -= vol
                    sample[well] = vol
                    vol -=vol
                    break
                else:
                    AvailableVol = self.wells[well]-self.MinWellVol
                    sample[well] = AvailableVol
                    self.wells[well] -= AvailableVol
                    vol -= AvailableVol
            else:
                pass # next well
        if vol !=0:
            print('{}: \tVol not reached'.format(self.compound))
        return sample

class Block():
    def __init__(self, Compound,DMSO, WorkingVol):
        self.WorkingVol = WorkingVol *1000 # convert to nl
        self.ProteinConc = 10 #ish
        self.K = 3 # prevents duplicates of zero values at 20 ul working vol
        self.Percent_DMSO = 0.05 # as a fraction of 1
        self.Compound = Compound
        self.DMSO = DMSO
        self.TestWells = ['X'+str(i) for i in range(1,9)]
        self.BlankWells = ['Y'+str(i) for i in range(1,9)]
        
        self.Transfers = self.MakeTransfer()
        
        
    def MakeTransfer(self):
        compound_vol = np.linspace(0,1,8)**self.K
        compound_vol *= self.Percent_DMSO
        compound_vol *= self.WorkingVol
        compound_vol = 2.5* np.round(compound_vol/2.5)
        TotalDMSOVol = np.round((self.Percent_DMSO * self.WorkingVol)/2.5) *2.5
        DMSO = TotalDMSOVol - compound_vol

        vols = {self.Compound:[i for i in compound_vol],\
             self.DMSO:[i for i in DMSO]}
        
        output = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        for vol_cpd, vol_dmso ,testwell, blankwell in zip(vols[self.Compound], vols[self.DMSO],self.TestWells, self.BlankWells):
            cpd_transfer = self.Compound.Sample(vol_cpd)
            dmso_transfer = self.DMSO.Sample(vol_dmso)
            for i,j in zip(cpd_transfer, dmso_transfer):
                temp = pd.DataFrame(\
                [[i,testwell,cpd_transfer[i]],\
                [j,testwell,dmso_transfer[j]],\
                [i,blankwell,cpd_transfer[i]],\
                [j,blankwell,dmso_transfer[j]]],\
                columns = ['SrcWell','Destination Well','Volume'])
                output = output.append(temp)
        return output.reset_index(drop=True)
            
class AssayPlate():
    def __init__(self):
        self.blocks = {}
        self.Allocations = {} # IDs map to blocks
        self.TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        
    def AddBlocks(self,block):
        count = len(self.blocks)+1
        self.blocks[count] = block

    def MapWells(self):
        # Gets the block from the number and replaces the
        # place holder destingation well IDs with actual well numbers
        TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        alphabet = string.ascii_uppercase
        for BlockNumber in self.blocks:
            # if odd
            if BlockNumber%2 == 1: 
                TestWells = [i+str(BlockNumber) for i in list(alphabet)[0:8]]
                BlankWells = [i+str(BlockNumber+1) for i in list(alphabet)[0:8]]
            else:
                TestWells = [i+str(BlockNumber-1) for i in list(alphabet)[8:16]]
                BlankWells = [i+str(BlockNumber) for i in list(alphabet)[8:16]]
                

            mappings = dict(zip(['X'+str(i) for i in range(1,9)]+['Y'+str(i) for i in range(1,9)],\
            TestWells + BlankWells))
            transfers = self.blocks[BlockNumber].Transfers
            transfers['Destination Well'] = transfers['Destination Well'].replace(mappings)
            self.TransferPlan = self.TransferPlan.append(transfers)
            self.TransferPlan.reset_index(inplace=True,drop=True)
    
    
   
def test_1():
    aracadonic = SourcePlateCompound('Arachadonic acid',['A'+str(i) for i in range(1,8)],ldv=False)
    DMSO = SourcePlateCompound('DMSO',['B'+str(i) for i in range(1,15)],ldv=False)
    assayplate = AssayPlate()
    for vol in [20,30,40]:
        for repeat in range(8):
            block = Block(aracadonic,DMSO,vol)
            assayplate.AddBlocks(block)
    assayplate.MapWells()
    transferMap = assayplate.TransferPlan
    transferMap = transferMap.loc[transferMap['Volume']>0]
    
    for i in transferMap['Destination Well'].unique():
        print(i, len(transferMap.loc[transferMap['Destination Well'] == i]))

def test_2():
    aracadonic = SourcePlateCompound('Arachadonic acid',['A'+str(i) for i in range(1,20)],ldv=False)
    DMSO = SourcePlateCompound('DMSO',['B'+str(i) for i in range(1,20)],ldv=False)

    assayplate1 = AssayPlate()
    assayplate2 = AssayPlate()
    assayplate3 = AssayPlate()
    assayplate4 = AssayPlate()
    
    for vol in [20,30,40]:
        for repeat in range(8):
            block1 = Block(aracadonic,DMSO,vol)
            assayplate1.AddBlocks(block1)
            block2 = Block(aracadonic,DMSO,vol)
            assayplate2.AddBlocks(block2)
            block3 = Block(aracadonic,DMSO,vol)
            assayplate3.AddBlocks(block3)
            block4 = Block(aracadonic,DMSO,vol)
            assayplate4.AddBlocks(block4)

    assayplate1.MapWells()
    assayplate2.MapWells()
    assayplate3.MapWells()
    assayplate4.MapWells()

    transferMap1 = assayplate1.TransferPlan
    transferMap2 = assayplate2.TransferPlan
    transferMap3 = assayplate3.TransferPlan
    transferMap4 = assayplate4.TransferPlan
    
    
if __name__ == '__main__':
    test_2()

import numpy as np
import pandas as pd
import plates

def main():
    df = pd.read_csv('compounds.csv') 
    smiles = df['SMILES'].to_list() + ['none', 'none']
    names = df['Item Name'].tolist() + ['none', 'none']


    concs = np.array([500 / 2**i for i in range(1,9)][::-1])
    
    blank = plates.UV384m3('data/blank.csv', 
                parser = 'ascii')
    wt = plates.UV384m3('data/3march-echo-wt-assay.csv', 
                parser = 'ascii')
    dm = plates.UV384m3('data/dm.csv', 
                parser = 'ascii')
    dm.report(save_dir = 'dm', 
                controlPlate = None, 
                concs = concs, 
                smiles = smiles, 
                names = names) 

    wt.report(save_dir = 'wt', 
                controlPlate = None, 
                concs = concs, 
                smiles = smiles, 
                names = names) 

if __name__ == '__main__':
    main()
import uv

data = uv.BM3('20210303-concs.csv')
print(data.concs)
concs = data.concs / 0.005

def v1(c1, c2, v2):
    return (c2 * v2) / c1 

print(v1(concs, 10, 30_000)) # µM, µM, µl
import numpy as np
import pandas as pd
import plates

def main():
    layout = pd.read_csv('compounds.csv', index_col=0)
    smiles = layout['SMILES']
    names = layout['Item Name']
    concs = np.array([500 / 2**i for i in range(1,9)])
    blank_data = plates.UV384m2('data/blank-pipette.CSV')
    wt_data = plates.UV384m2('data/wt-pipette.CSV') 
    dm_data = plates.UV384m2('data/dm-pipette.csv') 
    wt_data.report(smiles = smiles, names = names, concs = concs)
    dm_data.report(smiles = smiles, names = names, concs = concs)

if __name__ == '__main__':
    main()
import os
import numpy as np
import pandas as pd
import plates

def main():
    # smiles and names
    df = pd.read_csv('compounds.csv', index_col=0)
    d12_names = df['Item Name'][:24]
    d12_smiles = df['SMILES'][:24]
    d34_names = df['Item Name'][25:49]
    d34_smiles = df['SMILES'][25:49]
    concs = np.array([500 / 2**i for i in range(1,9)][::-1])

    # reports
    dfs = []
    for i, j in zip(['wt-d1.CSV','a82-d1,2.CSV','dm-d1,2.CSV'], ['wt','a82f','a82f-f87v']):
        plate = plates.UV384m1(i,
                name = j + 'd12', 
                concs = concs, 
                smiles = d12_smiles, 
                names = d12_names)
        dfs.append(plates.reportPlate(plate, save_dir = j + 'd12'))

    for i, j in zip(['wt-d3+4.CSV','a82-d3,4.CSV','dm-d3,4.CSV'], ['wt','a82f','a82f-f87v']):
        plate = plates.UV384m1(i,
                name = j + '34', 
                concs = concs, 
                smiles = d34_smiles, 
                names = d34_names)
        dfs.append(plates.reportPlate(plate, save_dir = j + 'd34'))

    df = pd.concat(dfs).reset_index(drop=True)
    df.to_csv('screeningSummary.csv')

    # anomaly detection
    # csv aggregation

if __name__ == '__main__':
    main()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class dataset:
    def __init__(self, path):
        self.data = pd.read_csv(path)
        self.headers = self.data.columns
        self.wavelengths = self.Get_Wavelengths()
        self.data = self.Get_numericalData() # throws away the metadata
        self.data = self.Zero_at_800()
        self.conc_labels = self.calc_conc()

    def Get_Wavelengths(self):
        wavelengths = self.data.iloc[:,0] # first column
        wavelengths = wavelengths[wavelengths.str.contains(r'\d\d\d.\d\d\d\d')].astype(float)
        # there's an integer in there somewhere!
        wavelengths = wavelengths.reset_index(drop=True).iloc[1:]
        # cba to sort this out now, I'm just going to plot it out of shift by one cell shit
        return wavelengths.reset_index(drop=True)

    def Get_numericalData(self):
        data = self.data
        data.columns = data.iloc[0,:]

        data = data.iloc[self.wavelengths.index,:].dropna(axis = 1)
        data = data.drop('Wavelength (nm)', axis = 1)
        data = data.iloc[1:,:] #drops strinds
        data.reset_index(inplace=True,drop=True)
        return data

    def Zero_at_800(self):
        data = self.data.astype(float)
        data = data.transpose()
        data.columns = self.wavelengths[:-1]
        zero_vals = data.iloc[:,0] # starts with 800
        data = data.subtract(zero_vals,axis=0)
        return data

    def calc_conc(self):
        data = self.data
        data = data.loc[:,data.columns < 421]
        data = data.loc[:,data.columns > 419]
        A420 = data.iloc[:,0]
        ext = 95
        conc_mM = A420/ext
        conc_uM = conc_mM*1000
        conc_uM.index = self.headers[::2][:-1]
        conc_uM.name = 'Concs'
        return conc_uM

    def plot_traces(self):
        data = self.data
        fig, ax = plt.subplots(figsize=(15,5))
        ax.set_prop_cycle('color',plt.cm.viridis(np.linspace(0,0.9,len(data))))
        for i in range(len(data)):
            y = data.iloc[i,:]
            plt.plot(y, lw = 4)
        plt.xlim((250,800))
        plt.ylim((-0.1,1))
        plt.xticks(np.linspace(250,800,11))
        plt.xlabel('Wavlength nm')
        plt.ylabel('Absorbance')
        plt.legend(self.conc_labels.round(2), title='P450 BM3 conc/uM')
        plt.title('P450 BM3 Wild Type Heme domain Concentration Check - 5th June 2019')
        plt.show()

def test():
    D = dataset('20191217_bm3WT_workingConcs.csv')
    print(D.calc_conc())
    
import sys 
import re 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from rdkit import Chem
from rdkit.Chem import Draw

class Plate:
    def __init__(self, csv):
        self.path = csv 
    @property
    def meta(self):
        # dict
        with open(self.path,'r') as f:
             m = ''.join([f.readline() for i in range(6)])
        metadata = {'date': re.findall('Date: (\d+/\d+/\d+)', m)[0],
                   'time': re.findall('Time: (\d+:\d+:\d+)', m)[0],
                   'protocol':re.findall('Test name: (.+?),',m)[0],
                   'plateID':re.findall('ID1: (.+)', m)[0]}
        return metadata

    @property
    def df(self):
        # clean
        df = pd.read_csv(self.path, skiprows=7)
        df.index = df.iloc[:,0] # wells
        df.drop(['Unnamed: 0','Wavelength [nm]'],axis=1, inplace=True)
        df.dropna(axis=1, inplace=True)
        df.columns = df.columns.astype(int) # wavelengths
        return df.astype(float)
    def row(self, letter):
        idx = pd.Series(self.df.index).str.contains(letter.upper())
        return self.df.loc[idx.to_list(), :]
    def col(self, num):
        idx = [int(re.findall('[A-Z]([0-9]+)', i)[0]) == int(num)  for i in self.df.index]
        return self.df.loc[idx, :]


class AssayPlate(Plate):
    def __init__(self, path):
        super().__init__(path)
    
    @property
    def blocks(self):
        return list(set([int(re.findall('[A-Z]([0-9]+)', i)[0])  for i in self.df.index]))
    def block(self, num):
        col = self.col(num)
        return col[::2], col[1::2] # samples, controls

    def normalizeBlock(self, num, smooth = False):
        # block number
        samples, controls = [i.subtract(i[800], axis = 0).reset_index(drop=True) for i in self.block(num)] # normalize at 800 nm
        samples = samples.subtract(controls, axis = 0)
        # scale at 405 nm inflection point
        inflection = samples.loc[:,405]
        scaling = 1 / inflection
        normalizedSamples = samples.multiply(scaling, axis = 0)
        if smooth:
            normalizedSamples = gaussian_filter1d(normalizedSamples, 1, axis = 1) # gaussian smoothing
            normalizedSamples = pd.DataFrame(normalizedSamples, index = samples.index, columns = samples.columns)
        return normalizedSamples

    def differenceBlock(self, num):
        # block num
        data = self.normalizeBlock(num)
        return data.subtract(data.loc[0,:], axis = 1)

    def responseBlock(self, num):
        # block num
        data = self.differenceBlock(num)
        return data.loc[:,420] - data.loc[:, 390]

def subplotTraces(data, ax, concs = None):
    # dataframe format: Plate.df
    # concs - array-like
    if concs is None:
        colors = plt.cm.inferno(np.linspace(0,1,len(data)))
    else:
        assert len(concs) == len(data)
        minmaxscale = lambda x : (x - min(x)) / max(x)
        colors = plt.cm.inferno(minmaxscale(concs))
    for i, j in zip(data.index, colors):
        ax.plot(data.loc[i, 300:], c = j, lw = 1)
    if concs is None:
        ax.legend(range(len(data)), title = 'Trace Number', loc = 'right')
    else:
        ax.legend(range(len(data)), labels = list(concs), title = '[Ligand] µM', loc = 'right')
    ax.set_xlabel('Wavelength nm')
    ax.set_ylabel('Absorbance')
    ax.set_ylim(-0.5,2)
    ax.set_xlim(300, 800)
    ax.set_xticks(np.linspace(300, 800, 11))
    ax.vlines(390,-0.5,2, linestyle = '--', lw = 0.5, color = '0.5')
    ax.vlines(420,-0.5,2, linestyle = '--', lw = 0.5, color = '0.5')

def subplotMichaelisMenten(x, y, ax):
    km, vmax = MichaelisMenten(x,y)
    ax.scatter(x,y)
    mm = lambda x, km, vmax : (x * vmax) / (km + x)
    xx = np.linspace(min(x), max(x), 100)
    ax.plot(xx, mm(xx, km, vmax))
    ax.set_ylim(min(y) - 0.1, max(y) * 1.1)
    ax.set_xlabel('[Ligand] µM')
    ax.set_ylabel('response')

def subplotText(ax, dictionary):
    s = ''.join([f'{i} = {round(j, 3)}\n' if type(j) == float else f'{i} = {j}\n' for i, j in zip(dictionary.keys(), dictionary.values())])
    ax.text(0.5,0.5, s, ha='center')
    ax.axis('off')


def MichaelisMenten(x,y):
    mm = lambda x, km, vmax : (x * vmax) / (km + x)
    (km, vmax), covariance = curve_fit(mm, x, y, bounds=((0,0),(np.inf, np.inf)))
    return km, vmax

def report(plate, num, concs = None, name = None, smiles = None, save_path = None):
    if concs is None:
        concs= np.array(range(8))
    fig = plt.figure(figsize = (10,10))
    if smiles is None:
        grid = plt.GridSpec(2,2)
    else:
        grid = plt.GridSpec(3,2)

    ax1 = fig.add_subplot(grid[0,0])
    subplotTraces(plate.normalizeBlock(num, smooth = True), ax1, concs)

    ax2 = fig.add_subplot(grid[1,0])
    subplotTraces(plate.differenceBlock(num), ax2, concs)

    ax3 = fig.add_subplot(grid[0,1])
    subplotMichaelisMenten(concs, plate.responseBlock(num), ax3)

    ax4 = fig.add_subplot(grid[1,1])
    km, vmax = MichaelisMenten(concs, plate.responseBlock(num))
    labels = {'km':round(km, 3) ,'vmax': round(vmax, 3)}
    if name is not None:
        labels['Name'] = name
    subplotText(ax4, labels)
    
    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        ax5 = fig.add_subplot(grid[2,0])
        im = Draw.MolToImage(mol)
        ax5.imshow(im)
        ax5.axis('off')
        
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)
    else:
        plt.show()
import pandas as pd
import matplotlib.pyplot as plt
import re

class UV:
    def __init__(self, path):
        self.path = path
        # slicing / dict indexing would be cool
    
    @property
    def df(self):
        return clean(self.path)

    def clean(self, path):
        df = pd.read_csv(path) # dataframe object from csv
        headers = df.columns # save column headers
        df = df.iloc[1:,:] # trim top row
        df.columns = headers # replace with old headers
        ### to do: #####
        ## make wavelengths index
        df.index = df.iloc[:,0]
        ## remove wavelength cols
        df = df.iloc[:,1::2]
        ## remove machine info (remove nan rows)
        df = df.dropna()
        ## get sample names from headers
        ## col headers to sample names
        df.columns = headers[::2][:-1]
        # round wavelenths
        df.index = [round(float(i)) for i in df.index]
        return df

    def plot_traces(self, df):
        plt.figure(figsize=(15,7)) # manually make canvas
        for col in df: # loop through columns
            plt.plot(df[col], # plot columns
                     label=col) # set trace label as col name - for legend

        plt.legend(loc='right') # detects 'label' in plot
        plt.xlabel('wavelength nm')
        plt.ylabel('absorbance')
        plt.ylim(-0.1,1)
        plt.xlim(250,800)
        plt.xticks(range(250,800,50)) # x axis ticks every 50 nm
        plt.title('UV-Vis absorbance of P450 BM3 with additions of arachadionic acid')
        plt.show()

    def regex_substrate_vols_ans(self, cols):
        # find only column headers with 'Arachadonic_acid'
        Arachadonic_acid_cols = [i for i in cols if 'acid' in i]
        # from those, extract numbers (volume of arachadionic acid added in µl)
        vols = [re.search('\d+\.\d',i).group() for i in Arachadonic_acid_cols]
        # also works
        # vols = [re.findall('\d+\.\d',i)[0] for i in Arachadonic_acid_cols]
        # convert the strings of the numebrs to floats
        vols = [float(i) for i in vols]
        # return the vol
        return vols

    def calc_concentrations_ans(self, v1,v2, c1):
        return (c1*v1)/v2

    def calc_change_a420_a390_ans(self, x):
        x = x.subtract(x.iloc[:,0], axis =0)
        a420 = x.loc[420,:]
        a390 = x.loc[390]
        total_change = a420.abs() + a390.abs()
        return total_change
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

import re

from sklearn.metrics import r2_score

from scipy import ndimage
from scipy.optimize import curve_fit

def argParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs=1, help='Trace')
    parser.add_argument('-s', nargs=1, help='Save')
    args=parser.parse_args()
    return args

class SpecData:
    def __init__(self, path,save=False):
        self.path=path
        self.save=save
        self.data = pd.read_csv(path)
        self.headers = pd.Series(self.data.columns[0::2])[:-1] # series to make it more workable #last col is blank
        self.substrateCols = self.regex_substrates()
        self.substrateConcs = self.regex_substrateConcs(self.substrateCols)

        self.wavelengths = self.Get_Wavelengths()
        self.data = self.Get_numericalData() # throws away the metadata
        self.data.index = self.headers
        self.data = self.Zero_at_800()
        self.data = self.data.iloc[self.substrateCols.index,:]
        self.data = self.GaussianFilterData()
        self.Diff = self.calcDiff(self.data)
        self.DiffDiff = self.calcDiffDiff()

    def Get_Wavelengths(self):
        wavelengths = self.data.iloc[:,0].dropna() # first column
        wavelengths = wavelengths[wavelengths.str.contains(r'\d\d\d.\d\d\d\d').dropna()].astype(float)
        # there's an integer in there somewhere!
        wavelengths = wavelengths.reset_index(drop=True).iloc[1:]
        # cba to sort this out now, I'm just going to plot it out of shift by one cell shit
        return wavelengths.reset_index(drop=True)

    def Get_numericalData(self):
        data = self.data
        data.columns = data.iloc[0,:]

        data = data.iloc[self.wavelengths.index,:].dropna(axis = 1)
        data = data.drop('Wavelength (nm)', axis = 1)
        data = data.iloc[1:,:] #drops strinds
        data.reset_index(inplace=True,drop=True)
        return data.transpose()

    def Zero_at_800(self):
        data = self.data.astype(float)
        data.columns = self.wavelengths[:-1]
        zero_vals = data.iloc[0,:] # starts with 800
        data = data.subtract(zero_vals,axis=1)
        return data

    def calc_conc(self):
        data = self.data
        data = data.loc[:,data.columns < 421]
        data = data.loc[:,data.columns > 419]
        A420 = data.iloc[:,0]
        ext = 95
        conc_mM = A420/ext
        conc_uM = conc_mM*1000
        conc_uM.reset_index(drop=True,inplace=True)
        conc_uM.name = 'P450 conc/uM'
        conc_uM.index=self.headers[:-1]
        return conc_uM

    def plot_traces(self, data,title=None):
        concs = self.substrateConcs
        fig, ax = plt.subplots(figsize=(10,6))
        ax.set_prop_cycle('color',plt.cm.inferno(np.linspace(0,0.9,len(data))))
        for i in range(len(data)):
            y = data.iloc[i,:]
            plt.plot(y, lw = 2)
        plt.xlim((250,800))
        plt.ylim((-0.1,0.6))
        plt.xticks(np.linspace(250,800,11))
        plt.xlabel('Wavlength nm')
        plt.ylabel('Absorbance')
        plt.legend(concs, title='Substrate conc/uM',loc = 'right')
        if title != None:
            plt.title(title)
        plt.show()

    def PlotMichaelesMenten(self,title = None):
        concs = self.substrateConcs
        km, vmax, loss = self.FitMichaelisMenten()
        km, vmax, loss  = km.item(), vmax.item(), loss.item()
        x = np.linspace(concs.min(),concs.max(),100)

        plt.figure(figsize=(7,7))
        plt.scatter(concs,self.DiffDiff)

        plt.plot(x, (vmax*x)/(km + x))
        plt.title(title)
        plt.xlabel('Concentration uM')
        plt.ylabel('Change in Absorbance')

        plt.text(concs.max()*0.7,vmax*0.2,'Km = '+str(np.around(km,2))+'\n'\
        +'Vmax = '+str(np.around(vmax,2))+'\n'\
        +'Loss = '+str(np.around(loss,5))+'\n'\
        'R squared = ' + str(round(1.-loss,3)))

        plt.show()


    def CalcRZ(self):
        #copy of the data
        data=self.data
        data.columns=self.wavelengths[:-1].round(0).astype(int)
        data.index=self.headers[:-1]
        RZ=self.data.loc[:,420]/self.data.loc[:,280]
        return RZ

    def regex_substrates(self):
        #substrates = self.headers.loc[self.headers.str.contains(expression, flags=re.IGNORECASE)]
        dmso = self.headers.loc[self.headers.str.contains('dmso', flags=re.IGNORECASE)]
        substrates = self.headers.loc[dmso.index.max()+1:]
        ProteinDMSO = pd.Series(self.headers.loc[(substrates.index.min() -1)])
        substrateCols = ProteinDMSO.append(substrates)

        # sort out index - DMSO needs to have index substrates.index.min() -1
        idx = list(substrateCols.index)
        idx[0] = substrates.index.min() -1
        substrateCols.index = idx
        return substrateCols

    def regex_substrateConcs(self,col_headers):
        vols =  col_headers[1:].str.extract('_(\d+(\.\d+)?)').astype(float)[0] # uls
        concs = (10_000*vols)/1000 # uM and ul
        concs = pd.Series([0.0]).append(concs)
        return concs

    def GaussianFilterData(self):
        data = pd.DataFrame(ndimage.gaussian_filter(self.data,2))
        data.columns = self.wavelengths[:-1]
        return data

    def calcDiff(self,data):
        return data-data.iloc[0,:]

    def calcDiffDiff(self):
        diff = self.Diff

        sec420 = diff.loc[:,diff.columns>416]
        sec420 = sec420.loc[:,sec420.columns<421]
        sec420 = sec420.loc[:,sec420.sum(axis=0).idxmax()]

        sec390 = diff.loc[:,diff.columns>385]
        sec390 = sec390.loc[:,sec390.columns<395]
        sec390 = sec390.loc[:,sec390.sum(axis=0).idxmax()]

        return (sec390 - sec420).abs()

    def FitMichaelisMenten(self):
        def michaelis_menten(x,km,vmax):
            return (vmax*x)/(km + x)
        params, covariance = curve_fit(michaelis_menten,
        self.substrateConcs,
        self.DiffDiff)
        km, vmax = params[0], params[1]
        pred = michaelis_menten(self.substrateConcs,km,vmax)
        r2 = r2_score(self.DiffDiff, pred)
        return km, vmax, r2

        '''
        x = self.substrateConcs.values
        x = torch.tensor(x,dtype=torch.float)

        y = self.DiffDiff.values
        y = torch.tensor(y,dtype=torch.float)

        km = torch.tensor([0.5],requires_grad=True,dtype=torch.float)
        vmax = torch.tensor([0.5],requires_grad=True,dtype=torch.float)
        optimizer = torch.optim.Adam({km,vmax},lr = 1e-2)
        loss_fn =  self.r_squared_torch

        for i in tqdm(range(5_000)):
            y_pred = (vmax*x)/(km + x) # scaling
            loss = 1 - loss_fn(y,y_pred) # 1 - r squared
            if km <0:
                # making sure that Km isn't negative
                loss -= km.item()
            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
        return km, vmax, loss'''

    def r_squared_torch(self,y,yh):
        residuals = y-yh
        ss_res = (residuals**2).sum()
        ss_tot = (y-y.mean()**2).sum()
        r_squared = 1 - (ss_res / ss_tot)
        return r_squared


def main():
    args=argParser()
    path=args.i[0]

    Dataset = SpecData(path)
    Dataset.plot_traces(Dataset.data, 'UV-Vis absorbance of Arachadonic Acid Titration into P450 BM3 WT')
    Dataset.plot_traces(Dataset.Diff, 'Relative UV-vis shift of Arachadonic Acid Titration into P450 BM3 WT')
    Dataset.PlotMichaelesMenten('Response of P450 BM3 WT to Arachadonic Acid Titration \n and Michaelis-Menten Curve Fit ')


if __name__ == '__main__':
    main()
import os
import numpy as np
import matplotlib.pyplot as plt
from plates import AssayPlate, subplotTraces, subplotMichaelisMenten, MichaelisMenten, subplotText, report
from tqdm import tqdm 


def main():
    os.makedirs('test-reports', exist_ok=True)
    plate = AssayPlate('a82-d1,2.CSV')
    for i in tqdm(plate.blocks):
        report(plate, i, name = '?', smiles = 'c1ccccc1', save_path=f'test-reports/block{i}.png')

if __name__ == '__main__':
    main()
import pandas as pd
import matplotlib.pyplot as plt
import re


def clean_data_ans(path):
    df = pd.read_csv(path) # dataframe object from csv
    headers = df.columns # save column headers
    df = df.iloc[1:,:] # trim top row
    df.columns = headers # replace with old headers
    ### to do: #####
    ## make wavelengths index
    df.index = df.iloc[:,0]
    ## remove wavelength cols
    df = df.iloc[:,1::2]
    ## remove machine info (remove nan rows)
    df = df.dropna()
    ## get sample names from headers
    ## col headers to sample names
    df.columns = headers[::2][:-1]
    # round wavelenths
    df.index = [round(float(i)) for i in df.index]
    return df

def plot_traces(df):
    plt.figure(figsize=(15,7)) # manually make canvas
    for col in df: # loop through columns
        plt.plot(df[col], # plot columns
                 label=col) # set trace label as col name - for legend

    plt.legend(loc='right') # detects 'label' in plot
    plt.xlabel('wavelength nm')
    plt.ylabel('absorbance')
    plt.ylim(-0.1,1)
    plt.xlim(250,800)
    plt.xticks(range(250,800,50)) # x axis ticks every 50 nm
    plt.title('UV-Vis absorbance of P450 BM3 with additions of arachadionic acid')
    plt.show()

def regex_substrate_vols_ans(cols):
    # find only column headers with 'Arachadonic_acid'
    Arachadonic_acid_cols = [i for i in cols if 'acid' in i]
    # from those, extract numbers (volume of arachadionic acid added in µl)
    vols = [re.search('\d+\.\d',i).group() for i in Arachadonic_acid_cols]
    # also works
    # vols = [re.findall('\d+\.\d',i)[0] for i in Arachadonic_acid_cols]
    # convert the strings of the numebrs to floats
    vols = [float(i) for i in vols]
    # return the vol
    return vols

def calc_concentrations_ans(v1,v2, c1):
    return (c1*v1)/v2

def calc_change_a420_a390_ans(x):
    x = x.subtract(x.iloc[:,0], axis =0)
    a420 = x.loc[420,:]
    a390 = x.loc[390]
    total_change = a420.abs() + a390.abs()
    return total_change
import pandas as pd
from picklist import SourcePlateCompound, Block, AssayPlate

def main():
    df = pd.read_csv('round1-cpds.csv')
    print(df)


if __name__ == '__main__':
    main()
import string
import numpy as np
import pandas as pd
from tqdm import tqdm

class SourcePlateCompound:
    def __init__(self, coumpound_name, wells, plate_id, well_vol = None, ldv = True):
        # give vol in µl
        assert isinstance(wells, list)
        self.compound = coumpound_name
        self.ldv = ldv
        self.plate_id = plate_id
        if self.ldv:
            self.MaxWellVol = 12 * 1000 #nl
            self.MinWellVol = 2.5 * 1000 #nl + safety
        else:
            self.MaxWellVol = 65 * 1000 #nl
            self.MinWellVol = 15 * 1000 #nl + safety

        if well_vol is None:
            self.well_vol = self.MaxWellVol
        else:
            self.well_vol = well_vol * 1000 # µl -> nl
        assert self.well_vol <= self.MaxWellVol
        self.wells = self.FillWells(wells)

    def FillWells(self,wells):
        output = {}
        for i in wells:
            output[i] = self.well_vol
        return output

    @property
    def AvailableVolume(self):
        return sum(self.wells.values())

    def Sample(self,vol):
        if vol %2.5 !=0:
            print('Transfer vol not a multiple of 2.5')
        sample = {}
        # take what you can then move on to the next well
        for well in self.wells:
            if self.wells[well] > self.MinWellVol:
                if vol < self.wells[well] - self.MinWellVol:
                    self.wells[well] -= vol
                    sample[well] = vol
                    vol -=vol
                    break
                else:
                    AvailableVol = self.wells[well]-self.MinWellVol
                    sample[well] = AvailableVol
                    self.wells[well] -= AvailableVol
                    vol -= AvailableVol
            else:
                pass # next well
        if vol !=0:
            print(f'{self.compound}: \tVol not reached \t {self.AvailableVolume / 1000} µl')
        return sample

class Block():
    def __init__(self, Compound, DMSO, WorkingVol):
        self.WorkingVol = WorkingVol *1000 # convert to nl
        self.ProteinConc = 10 #ish
        self.K = 3 # prevents duplicates of zero values at 20 ul working vol
        self.Percent_DMSO = 0.05 # as a fraction of 1
        self.Compound = Compound
        self.DMSO = DMSO
        self.TestWells = ['X'+str(i) for i in range(1,9)]
        self.Transfers = self.MakeTransfer()
        
    def MakeTransfer(self):
        compound_vol = np.linspace(0,1,8)**self.K
        compound_vol *= self.Percent_DMSO
        compound_vol *= self.WorkingVol
        compound_vol = 2.5* np.round(compound_vol/2.5)
        TotalDMSOVol = np.round((self.Percent_DMSO * self.WorkingVol)/2.5) *2.5
        DMSO = TotalDMSOVol - compound_vol

        vols = {self.Compound:[i for i in compound_vol],\
             self.DMSO:[i for i in DMSO]}
        # todo add source plate id 
        #   - is that a valid feild for echo picklist?
        output = pd.DataFrame([],columns = ['Plate',
                                        'SrcWell',
                                        'Destination Well',
                                        'Volume'])
        for vol_cpd, vol_dmso ,testwell in zip(vols[self.Compound], 
                                            vols[self.DMSO],
                                            self.TestWells):
            cpd_transfer = self.Compound.Sample(vol_cpd)
            dmso_transfer = self.DMSO.Sample(vol_dmso)
            # todo register source plate
            for i,j in zip(cpd_transfer, dmso_transfer):
                temp = pd.DataFrame(\
                [[self.Compound.plate_id, i,testwell,cpd_transfer[i]],\
                [self.DMSO.plate_id, j,testwell,dmso_transfer[j]]],\
                columns = ['Plate','SrcWell','Destination Well','Volume'])
                output = output.append(temp)
        # todo - split picklist or see if echo supports multiple input plates
        return output.reset_index(drop=True)
            
class AssayPlate():
    def __init__(self):
        self.blocks = {}
        self.TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        self.wells = [f'{i}{j}' for i in string.ascii_uppercase[:16] for j in range(1, 25)]
    def AddBlocks(self,block):
        count = len(self.blocks)+1
        self.blocks[count] = block

    def MapWells(self):
        TransferPlan = pd.DataFrame([],columns = ['Plate', 
                                    'SrcWell',
                                    'Destination Well',
                                    'Volume'])
        alphabet = string.ascii_uppercase
        for BlockNumber, wells in zip(self.blocks, 
                [self.wells[i*8:i*8 + 8] for i in range(round(len(self.wells)/8))]):
            mappings = dict(zip(['X'+str(i) for i in range(1,9)], wells))
            transfers = self.blocks[BlockNumber].Transfers
            transfers['Destination Well'] = transfers['Destination Well'].replace(mappings)
            self.TransferPlan = self.TransferPlan.append(transfers)
            self.TransferPlan.reset_index(inplace=True,drop=True)
        self.TransferPlan = self.TransferPlan.loc[self.TransferPlan['Volume'] > 0]
    
def main():
    df = pd.read_csv('round1-cpds.csv')

    dmso_sourcewells = [f'{i}{j}' for i in string.ascii_uppercase[:4] for j in range(1, 25)]
    DMSO = SourcePlateCompound('DMSO',dmso_sourcewells, 'A', well_vol = 60, ldv=False) # palte A

    cpd_src_wells = [f'{i}{j}' for i in string.ascii_uppercase[:16] for j in range(1, 25)]
    split = lambda l, n : [l[i*n:i*n + n] for i in range(round(len(l) / n))]
    cpd_src_tuple = split(cpd_src_wells, 4)
    CPDS = [SourcePlateCompound(i, j, 'B', well_vol = 11, ldv = True) for i, j in zip(df['Item Name'], cpd_src_tuple)]

    src = pd.DataFrame([[df.loc[i,'Item Name'], 
                        df.loc[i,'Rack Number'],
                        df.loc[i,'Plate Location'],
                        df.loc[i, 'CatalogNumber'],
                        k] for i, j in zip(df.index, cpd_src_tuple) for k in j],
                        columns = ['Item Name', 'Rack Number', 'Plate Location','CatalogNumber', 'Source Well'])
    src.to_csv('src-layout.csv')

    for repeat in tqdm(range(1,5)):
        assayplate1 = AssayPlate()
        for i in CPDS[:48]:
            block = Block(i,DMSO,30)
            assayplate1.AddBlocks(block)
        assayplate1.MapWells()
        transferMap1 = assayplate1.TransferPlan
        transferMap1.to_csv(f'A{repeat}.csv')

        assayplate2 = AssayPlate()
        for i in CPDS[:48]:
            block = Block(i,DMSO,30)
            assayplate2.AddBlocks(block)
        assayplate2.MapWells()
        transferMap2 = assayplate2.TransferPlan
        transferMap2.to_csv(f'B{repeat}.csv')

if __name__ == '__main__':
    main()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import string



class SourcePlateCompound():
    def __init__(self,coumpound_name,wells,ldv = True):
        self.compound = coumpound_name
        self.ldv = ldv
        if self.ldv:
            self.MaxWellVol = 12 * 1000 #nl
            self.MinWellVol = 2.5 * 1000 #nl + safety
        else:
            self.MaxWellVol = 65 * 1000 #nl
            self.MinWellVol = 20 * 1000 #nl + safety
        self.wells = self.FillWells(wells)

    def FillWells(self,wells,):
        output = {}
        
        if self.ldv:
            for i in wells:
                output[i] = 11 * 1000 #self.MaxWellVol #ul
        else:
            for i in wells:
                output[i] = 60 * 1000
        return output

    def AvailableVolume(self):
        return sum(self.wells.values())

    def Sample(self,vol):
        if vol %2.5 !=0:
            print('Transfer vol not a multiple of 2.5')
        sample = {}
        # take what you can then move on to the next well
        for well in self.wells:
            if self.wells[well] > (self.MinWellVol + 5000): # extra safety margin
                
                if vol < self.wells[well]:
                    self.wells[well] -= vol
                    sample[well] = vol
                    vol -=vol
                    break
                else:
                    AvailableVol = self.wells[well]-self.MinWellVol
                    sample[well] = AvailableVol
                    self.wells[well] -= AvailableVol
                    vol -= AvailableVol
            else:
                pass # next well
        if vol !=0:
            print('{}: \tVol not reached'.format(self.compound))
        return sample

class Block():
    def __init__(self, Compound,DMSO, WorkingVol):
        self.WorkingVol = WorkingVol *1000 # convert to nl
        self.ProteinConc = 10 #ish
        self.K = 3 # prevents duplicates of zero values at 20 ul working vol
        self.Percent_DMSO = 0.05 # as a fraction of 1
        self.Compound = Compound
        self.DMSO = DMSO
        self.TestWells = ['X'+str(i) for i in range(1,9)]
        self.BlankWells = ['Y'+str(i) for i in range(1,9)]
        
        self.Transfers = self.MakeTransfer()
        
        
    def MakeTransfer(self):
        compound_vol = np.linspace(0,1,8)**self.K
        compound_vol *= self.Percent_DMSO
        compound_vol *= self.WorkingVol
        compound_vol = 2.5* np.round(compound_vol/2.5)
        TotalDMSOVol = np.round((self.Percent_DMSO * self.WorkingVol)/2.5) *2.5
        DMSO = TotalDMSOVol - compound_vol

        vols = {self.Compound:[i for i in compound_vol],\
             self.DMSO:[i for i in DMSO]}
        
        output = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        for vol_cpd, vol_dmso ,testwell, blankwell in zip(vols[self.Compound], vols[self.DMSO],self.TestWells, self.BlankWells):
            cpd_transfer = self.Compound.Sample(vol_cpd)
            dmso_transfer = self.DMSO.Sample(vol_dmso)
            for i,j in zip(cpd_transfer, dmso_transfer):
                temp = pd.DataFrame(\
                [[i,testwell,cpd_transfer[i]],\
                [j,testwell,dmso_transfer[j]],\
                [i,blankwell,cpd_transfer[i]],\
                [j,blankwell,dmso_transfer[j]]],\
                columns = ['SrcWell','Destination Well','Volume'])
                output = output.append(temp)
        return output.reset_index(drop=True)
            
class AssayPlate():
    def __init__(self):
        self.blocks = {}
        self.Allocations = {} # IDs map to blocks
        self.TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        
    def AddBlocks(self,block):
        count = len(self.blocks)+1
        self.blocks[count] = block

    def MapWells(self):
        # Gets the block from the number and replaces the
        # place holder destingation well IDs with actual well numbers
        TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        alphabet = string.ascii_uppercase
        for BlockNumber in self.blocks:
            # if odd
            if BlockNumber%2 == 1: 
                TestWells = [i+str(BlockNumber) for i in list(alphabet)[0:8]]
                BlankWells = [i+str(BlockNumber+1) for i in list(alphabet)[0:8]]
            else:
                TestWells = [i+str(BlockNumber-1) for i in list(alphabet)[8:16]]
                BlankWells = [i+str(BlockNumber) for i in list(alphabet)[8:16]]
                

            mappings = dict(zip(['X'+str(i) for i in range(1,9)]+['Y'+str(i) for i in range(1,9)],\
            TestWells + BlankWells))
            transfers = self.blocks[BlockNumber].Transfers
            transfers['Destination Well'] = transfers['Destination Well'].replace(mappings)
            self.TransferPlan = self.TransferPlan.append(transfers)
            self.TransferPlan.reset_index(inplace=True,drop=True)
    
    
   
def test_1():
    aracadonic = SourcePlateCompound('Arachadonic acid',['A'+str(i) for i in range(1,8)],ldv=False)
    DMSO = SourcePlateCompound('DMSO',['B'+str(i) for i in range(1,15)],ldv=False)
    assayplate = AssayPlate()
    for vol in [20,30,40]:
        for repeat in range(8):
            block = Block(aracadonic,DMSO,vol)
            assayplate.AddBlocks(block)
    assayplate.MapWells()
    transferMap = assayplate.TransferPlan
    transferMap = transferMap.loc[transferMap['Volume']>0]
    
    for i in transferMap['Destination Well'].unique():
        print(i, len(transferMap.loc[transferMap['Destination Well'] == i]))

def test_2():
    aracadonic = SourcePlateCompound('Arachadonic acid',['A'+str(i) for i in range(1,20)],ldv=False)
    DMSO = SourcePlateCompound('DMSO',['B'+str(i) for i in range(1,20)],ldv=False)

    assayplate1 = AssayPlate()
    assayplate2 = AssayPlate()
    assayplate3 = AssayPlate()
    assayplate4 = AssayPlate()
    
    for vol in [20,30,40]:
        for repeat in range(8):
            block1 = Block(aracadonic,DMSO,vol)
            assayplate1.AddBlocks(block1)
            block2 = Block(aracadonic,DMSO,vol)
            assayplate2.AddBlocks(block2)
            block3 = Block(aracadonic,DMSO,vol)
            assayplate3.AddBlocks(block3)
            block4 = Block(aracadonic,DMSO,vol)
            assayplate4.AddBlocks(block4)

    assayplate1.MapWells()
    assayplate2.MapWells()
    assayplate3.MapWells()
    assayplate4.MapWells()

    transferMap1 = assayplate1.TransferPlan
    transferMap2 = assayplate2.TransferPlan
    transferMap3 = assayplate3.TransferPlan
    transferMap4 = assayplate4.TransferPlan
    
    
if __name__ == '__main__':
    test_2()

from string import ascii_uppercase
import pandas as pd
import echo

def main():
    wells = [f'{i}{j}' for i in ascii_uppercase[:16] for  j in range(1,25)]
    df = pd.read_csv('round1-cpds.csv')

    # dmso
    dmso = echo.Cpd(name = 'dmso', vol = 10_000)
    dmso_src = echo.Src(ldv = False)
    for i in wells[:24]:
        dmso_src.fill(i, dmso, 60)
    
    # compounds
    cpds = [echo.Cpd(name = i, vol = 20) for i in df['Item Name']]
    cpd_src = echo.Src()
    for i, j in zip(wells, cpds):
        cpd_src.fill(i, j)

    # blocks and dest fill
    dest = echo.Dest()
    # block1 = echo.Block(cpds[0], dmso, dest.wells) # slice?


if __name__ == '__main__':
    main()
import numpy as np
import pandas as pd
import plates

def main():
    concs = np.linspace(0,1,8)**3 * 500 # is k = 3?
    blank = plates.UV384m4('../data/blank-a.CSV')
    plate = plates.UV384m4('../data/a83f-a.CSV', 
            control = blank,
            concs = concs)
    for i in [1,2,3,4]:
        blk = plate.block(i, smiles = 'CCCCCC=O')
        #print(blk.df)
        #print(blk.norm)
        #print(blk.diff)
        #print(blk.response)
        print(blk.mm)
        #print(pd.DataFrame([plate.block(i).mm for i in range(1,24)]))
        #plates.report(blk)

if __name__ == '__main__':
    main()
import argparse
from tqdm import tqdm
import swalign

def parse(path):
    # assumes one sequence
    with open(path,'r') as f:
        data=f.read().splitlines()
    header = data[0]
    seq = ''.join(data[1:])
    return header, seq

def main(args):

    match=2
    mismatch=-1
    scoring=swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    
    wt_head, wt_seq = parse(args.ref)
    wt_seq = wt_seq.upper()

    with open(args.out,'w') as f:
        for i in tqdm(args.sequences):
            head, seq = parse(i)
            aln = sw.align(wt_seq, seq)
            f.write(f'{head} and {wt_head} \n')
            aln.dump(wrap = 60, out=f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--out', default='out')
    parser.add_argument('-r','--ref')
    parser.add_argument('-s','--sequences', nargs='+')
    args = parser.parse_args()
    main(args)
bm3= "MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQALKFVRDFAGDGLFTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASGEQSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPLGGIPSPSTEQSAKKVRK*"


orf = 'ATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCcatatgacaattaaagaaatgcctcagccaaaaacgtttggagagcttaaaaatttaccgttattaaacacagataaaccggttcaagctttgatgaaaattgcggatgaattaggagaaatctttaaattcgaggcgcctggccgtgtaacgcgctacttatcaagtcagcgtctaattaaagaagcatgcgatgaatcacgctttgataaaaacttaagtcaagcgcttaaatttgtacgtgattttgcaggagacgggttatttacaagctggacgcatgaaaaaaattggaaaaaagcgcataatatcttacttccaagcttcagtcagcaggcaatgaaaggctatcatgcgatgatggtcgatatcgccgtgcagcttgttcaaaagtgggagcgtctaaatgcagatgagcatattgaggtaccggaagacatgacacgtttaacgcttgatacaattggtctttgcggctttaactatcgctttaacagcttttaccgagatcagcctcatccatttattacaagtatggtccgtgcactggatgaagcaatgaacaagctgcagcgagcaaatccagacgacccagcttatgatgaaaacaagcgccagtttcaagaagatatcaaggtgatgaacgacctagtagataaaattattgcagatcgcaaagcaagcggtgaacaaagcgatgatttattaacgcacatgctaaacggaaaagatccagaaacgggtgagccgcttgatgacgagaacattcgctatcaaattattacattcttaattgcgggacacgaaacaactagtggtcttttatcatttgcgctgtatttcttagtgaaaaatccacatgtattacaaaaagcagcagaagaagcagcacgagttctagtagatcctgttccaagctacaaacaagtcaaacagcttaaatatgtcggcatggtcttaaacgaagcgctgcgcttatggccaactgctcctgcgttttccctatatgcaaaagaagatacggtgcttggaggagaatatcctttagaaaaaggcgacgaactaatggttctgattcctcagcttcaccgtgataaaacaatttggggagacgatgtggaagagttccgtccagagcgttttgaaaatccaagtgcgattccgcagcatgcgtttaaaccgtttggaaacggtcagcgtgcgtgtatcggtcagcagttcgctcttcatgaagcaacgctggtcctaggtatgatgctaaaacactttgactttgaagatcatacaaactacgagctggatattaaagaaactttaacgttaaaacctgaaggctttgtggtaaaagcaaaatcgaaaaaaattccgcttggcggtattccttcacctagcactgaacagtctgctaaaaaagtacgcaaatag'
import re
from tqdm import tqdm 
import pandas as pd 

import mxn
from bm3 import bm3, orf


COMPLIMENT = {'A':'T','T':'A','C':'G','G':'C'}

reverse_comp = lambda s : ''.join([COMPLIMENT[i.upper()] for i in s[::-1]])

def parse_mutation(m):
    return int(re.findall('\d+', m)[0]), re.findall('\w', m)[-1] 

def main():
    with open('mutations.txt', 'r') as f:
        mutations = f.read().split('\n')[:-1]

    output = []
    for i in tqdm(mutations):
        cds = mxn.CDS(orf, bm3)
        try:
            pos, aa = parse_mutation(i)
            cds.mutate(pos, aa)
            output.append(pd.DataFrame(cds.primers).T)
        except:
            print(i)
    df = pd.concat(output)
    df.to_csv('primers.csv')





if __name__ == '__main__':
    main()
import re
from tqdm import tqdm 
import numpy as np
import pandas as pd 
from scipy import optimize

import mxn
from mxn.agilent_tm import agilent_tm
from bm3 import bm3, orf

import primer3

COMPLIMENT = {'A':'T','T':'A','C':'G','G':'C'}

reverse_comp = lambda s : ''.join([COMPLIMENT[i.upper()] for i in s[::-1]])


def agilent_tm(s):
    mv_conc, dv_conc, dntp_conc, dna_conc, dmsoPerc = [25, # mv - 25-100
            0.5, # dv - 0.5 - 5
            0.8, # dntp
            1, # dna
            8] # dmso 0-10%
    tm = primer3.bindings.calcTm(s, 
            mv_conc = mv_conc, 
            dv_conc = dv_conc, 
            dntp_conc = dntp_conc, 
            dna_conc = dna_conc,
            tm_method = 'santalucia',
            salt_corrections_method = 'santalucia')
    return tm - 0.75 * dmsoPerc

def generic_primer(seq, comp='rev', tm=78):
    # seq = sense, region to make primer in
    if comp == 'rev':
        seq = reverse_comp(seq)

    homoTm = lambda s : primer3.calcHomodimerTm(s,mv_conc= 25, dv_conc = 0.5)
    endscore = lambda s : sum([1 for i in s[1] + s[-2] if i == 'C' or i == 'G']\
                            + [2 for i in s[0] + s[-1] if i == 'C' or i == 'G']) # max 6
    scorefn = lambda s : abs(agilent_tm(s) - tm) \
                        - (endscore(s) / 6) \
                        + max([(homoTm(s) / 6),0])
    select = lambda n1, n2 : seq[round(n1):round(n1)+round(n2)]
    objective = lambda n1, n2 : scorefn(select(n1, n2))
    helper = lambda array : objective(array[0], array[1])
    
    results = optimize.dual_annealing(helper, 
                            bounds = ((0,len(seq) - 60),
                                    (10,60)),
                            initial_temp = 1e5)
    primer = select(results.x[0], results.x[1])
    score = scorefn(primer)
    tm = agilent_tm(primer)
    return {'primer':primer,
            'tm':tm,
            'end_score': endscore(primer),
            'homotm':homoTm(primer),
            'length':len(primer)}

def main():
    primer = generic_primer(orf[-100:])





if __name__ == '__main__':
    main()
import numpy as np
import pandas as pd
import echo


def square_layout():
    # 
    rows = [echo.hwells[i:i+24] for i in range(0,384,24)]
    x=[]
    for i,j in zip(rows[::2],rows[1::2]):
        for k in range(0,24,2):
            x.append(i[k:k+2] + j[k:k+2])
    return x

def main():
    # plate 1 a1:a1
    # plate 2 a1:b1
    # plate 3 a1:a2
    # plate 4 a1:b2

    compounds = pd.read_csv('fda.csv', index_col=0)

    cpds = [echo.Cpd(i, vol=100) for i in compounds['CatalogNumber']]
    cpds_blocks_of_384 = [cpds[i:i+384] for i in range(0, len(cpds), 384)]
    dmso = echo.Cpd('dmso', vol=1500)

    src_dmso = echo.SrcPlate(name='src_dmso', ldv=False)


    src_plates = []
    for i,j in enumerate(cpds_blocks_of_384,1):
        src = echo.SrcPlate(name=f'src{i}')
        for k,l in zip(j, square_layout()):
            for m in l:
                src[m].fill(k.sample(12))
        src_plates.append(src)
    

    for i in src_dmso[192:192+19]:
        i.fill(dmso.sample(60))

    for plate in [1,2,3,4]:
        dest = echo.DestPlate(name=f'dest-{plate}')
        for i,j in zip(dest, cpds):
            j.xfer(12,i)


    #src_xfers = pd.DataFrame(src1.xfer_record)
    #dmso_xfers = pd.DataFrame(src_dmso.xfer_record)
    #src_plate_fill = pd.DataFrame(src_plate_fill)
    #src_plate_fill.to_csv('src_plate_fill.csv',index=False)
    #src_xfers.to_csv('cpd-src-picklist.csv')
    #dmso_xfers.to_csv('dmso-src-picklist.csv')


if __name__ == '__main__':
    main()
import pandas as pd
from picklist import SourcePlateCompound, Block, AssayPlate

def main():
    df = pd.read_csv('round1-cpds.csv')
    print(df)


if __name__ == '__main__':
    main()
import string
import numpy as np
import pandas as pd
from tqdm import tqdm

class SourcePlateCompound:
    def __init__(self, coumpound_name, wells, plate_id, well_vol = None, ldv = True):
        # give vol in µl
        assert isinstance(wells, list)
        self.compound = coumpound_name
        self.ldv = ldv
        self.plate_id = plate_id
        if self.ldv:
            self.MaxWellVol = 12 * 1000 #nl
            self.MinWellVol = 2.5 * 1000 #nl + safety
        else:
            self.MaxWellVol = 65 * 1000 #nl
            self.MinWellVol = 15 * 1000 #nl + safety

        if well_vol is None:
            self.well_vol = self.MaxWellVol
        else:
            self.well_vol = well_vol * 1000 # µl -> nl
        assert self.well_vol <= self.MaxWellVol
        self.wells = self.FillWells(wells)

    def FillWells(self,wells):
        output = {}
        for i in wells:
            output[i] = self.well_vol
        return output

    @property
    def AvailableVolume(self):
        return sum(self.wells.values())

    def Sample(self,vol):
        if vol %2.5 !=0:
            print('Transfer vol not a multiple of 2.5')
        sample = {}
        # take what you can then move on to the next well
        for well in self.wells:
            if self.wells[well] > self.MinWellVol:
                if vol < self.wells[well] - self.MinWellVol:
                    self.wells[well] -= vol
                    sample[well] = vol
                    vol -=vol
                    break
                else:
                    AvailableVol = self.wells[well]-self.MinWellVol
                    sample[well] = AvailableVol
                    self.wells[well] -= AvailableVol
                    vol -= AvailableVol
            else:
                pass # next well
        if vol !=0:
            print(f'{self.compound}: \tVol not reached \t {self.AvailableVolume / 1000} µl')
        return sample

class Block():
    def __init__(self, Compound, DMSO, WorkingVol):
        self.WorkingVol = WorkingVol *1000 # convert to nl
        self.ProteinConc = 10 #ish
        self.K = 3 # prevents duplicates of zero values at 20 ul working vol
        self.Percent_DMSO = 0.05 # as a fraction of 1
        self.Compound = Compound
        self.DMSO = DMSO
        self.TestWells = ['X'+str(i) for i in range(1,9)]
        self.Transfers = self.MakeTransfer()
        
    def MakeTransfer(self):
        compound_vol = np.linspace(0,1,8)**self.K
        compound_vol *= self.Percent_DMSO
        compound_vol *= self.WorkingVol
        compound_vol = 2.5* np.round(compound_vol/2.5)
        TotalDMSOVol = np.round((self.Percent_DMSO * self.WorkingVol)/2.5) *2.5
        DMSO = TotalDMSOVol - compound_vol

        vols = {self.Compound:[i for i in compound_vol],\
             self.DMSO:[i for i in DMSO]}
        # todo add source plate id 
        #   - is that a valid feild for echo picklist?
        output = pd.DataFrame([],columns = ['Plate',
                                        'SrcWell',
                                        'Destination Well',
                                        'Volume'])
        for vol_cpd, vol_dmso ,testwell in zip(vols[self.Compound], 
                                            vols[self.DMSO],
                                            self.TestWells):
            cpd_transfer = self.Compound.Sample(vol_cpd)
            dmso_transfer = self.DMSO.Sample(vol_dmso)
            # todo register source plate
            for i,j in zip(cpd_transfer, dmso_transfer):
                temp = pd.DataFrame(\
                [[self.Compound.plate_id, i,testwell,cpd_transfer[i]],\
                [self.DMSO.plate_id, j,testwell,dmso_transfer[j]]],\
                columns = ['Plate','SrcWell','Destination Well','Volume'])
                output = output.append(temp)
        # todo - split picklist or see if echo supports multiple input plates
        return output.reset_index(drop=True)
            
class AssayPlate():
    def __init__(self):
        self.blocks = {}
        self.TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        self.wells = [f'{i}{j}' for i in string.ascii_uppercase[:16] for j in range(1, 25)]
    def AddBlocks(self,block):
        count = len(self.blocks)+1
        self.blocks[count] = block

    def MapWells(self):
        TransferPlan = pd.DataFrame([],columns = ['Plate', 
                                    'SrcWell',
                                    'Destination Well',
                                    'Volume'])
        alphabet = string.ascii_uppercase
        for BlockNumber, wells in zip(self.blocks, 
                [self.wells[i*8:i*8 + 8] for i in range(round(len(self.wells)/8))]):
            mappings = dict(zip(['X'+str(i) for i in range(1,9)], wells))
            transfers = self.blocks[BlockNumber].Transfers
            transfers['Destination Well'] = transfers['Destination Well'].replace(mappings)
            self.TransferPlan = self.TransferPlan.append(transfers)
            self.TransferPlan.reset_index(inplace=True,drop=True)
        self.TransferPlan = self.TransferPlan.loc[self.TransferPlan['Volume'] > 0]
    
def main():
    df = pd.read_csv('round1-cpds.csv')

    dmso_sourcewells = [f'{i}{j}' for i in string.ascii_uppercase[5:9] for j in range(1, 25)]
    DMSO = SourcePlateCompound('DMSO',dmso_sourcewells, 'A', well_vol = 60, ldv=False) # palte A

    cpd_src_wells = [f'{i}{j}' for i in string.ascii_uppercase[:16] for j in range(1, 25)]
    split = lambda l, n : [l[i*n:i*n + n] for i in range(round(len(l) / n))]
    cpd_src_tuple = split(cpd_src_wells, 4)
    CPDS = [SourcePlateCompound(i, j, 'B', well_vol = 11, ldv = True) for i, j in zip(df['Item Name'], cpd_src_tuple)]

    src = pd.DataFrame([[df.loc[i,'Item Name'], 
                        df.loc[i,'Rack Number'],
                        df.loc[i,'Plate Location'],
                        df.loc[i, 'CatalogNumber'],
                        k] for i, j in zip(df.index, cpd_src_tuple) for k in j],
                        columns = ['Item Name', 'Rack Number', 'Plate Location','CatalogNumber', 'Source Well'])
    src.to_csv('src-layout.csv')

    for repeat in tqdm(range(1,5)):
        assayplate1 = AssayPlate()
        for i in CPDS[:48]:
            block = Block(i,DMSO,30)
            assayplate1.AddBlocks(block)
        assayplate1.MapWells()
        transferMap1 = assayplate1.TransferPlan
        transferMap1.to_csv(f'A{repeat}.csv')

        assayplate2 = AssayPlate()
        for i in CPDS[:48]:
            block = Block(i,DMSO,30)
            assayplate2.AddBlocks(block)
        assayplate2.MapWells()
        transferMap2 = assayplate2.TransferPlan
        transferMap2.to_csv(f'B{repeat}.csv')

if __name__ == '__main__':
    main()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import string



class SourcePlateCompound():
    def __init__(self,coumpound_name,wells,ldv = True):
        self.compound = coumpound_name
        self.ldv = ldv
        if self.ldv:
            self.MaxWellVol = 12 * 1000 #nl
            self.MinWellVol = 2.5 * 1000 #nl + safety
        else:
            self.MaxWellVol = 65 * 1000 #nl
            self.MinWellVol = 20 * 1000 #nl + safety
        self.wells = self.FillWells(wells)

    def FillWells(self,wells,):
        output = {}
        
        if self.ldv:
            for i in wells:
                output[i] = 11 * 1000 #self.MaxWellVol #ul
        else:
            for i in wells:
                output[i] = 60 * 1000
        return output

    def AvailableVolume(self):
        return sum(self.wells.values())

    def Sample(self,vol):
        if vol %2.5 !=0:
            print('Transfer vol not a multiple of 2.5')
        sample = {}
        # take what you can then move on to the next well
        for well in self.wells:
            if self.wells[well] > (self.MinWellVol + 5000): # extra safety margin
                
                if vol < self.wells[well]:
                    self.wells[well] -= vol
                    sample[well] = vol
                    vol -=vol
                    break
                else:
                    AvailableVol = self.wells[well]-self.MinWellVol
                    sample[well] = AvailableVol
                    self.wells[well] -= AvailableVol
                    vol -= AvailableVol
            else:
                pass # next well
        if vol !=0:
            print('{}: \tVol not reached'.format(self.compound))
        return sample

class Block():
    def __init__(self, Compound,DMSO, WorkingVol):
        self.WorkingVol = WorkingVol *1000 # convert to nl
        self.ProteinConc = 10 #ish
        self.K = 3 # prevents duplicates of zero values at 20 ul working vol
        self.Percent_DMSO = 0.05 # as a fraction of 1
        self.Compound = Compound
        self.DMSO = DMSO
        self.TestWells = ['X'+str(i) for i in range(1,9)]
        self.BlankWells = ['Y'+str(i) for i in range(1,9)]
        
        self.Transfers = self.MakeTransfer()
        
        
    def MakeTransfer(self):
        compound_vol = np.linspace(0,1,8)**self.K
        compound_vol *= self.Percent_DMSO
        compound_vol *= self.WorkingVol
        compound_vol = 2.5* np.round(compound_vol/2.5)
        TotalDMSOVol = np.round((self.Percent_DMSO * self.WorkingVol)/2.5) *2.5
        DMSO = TotalDMSOVol - compound_vol

        vols = {self.Compound:[i for i in compound_vol],\
             self.DMSO:[i for i in DMSO]}
        
        output = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        for vol_cpd, vol_dmso ,testwell, blankwell in zip(vols[self.Compound], vols[self.DMSO],self.TestWells, self.BlankWells):
            cpd_transfer = self.Compound.Sample(vol_cpd)
            dmso_transfer = self.DMSO.Sample(vol_dmso)
            for i,j in zip(cpd_transfer, dmso_transfer):
                temp = pd.DataFrame(\
                [[i,testwell,cpd_transfer[i]],\
                [j,testwell,dmso_transfer[j]],\
                [i,blankwell,cpd_transfer[i]],\
                [j,blankwell,dmso_transfer[j]]],\
                columns = ['SrcWell','Destination Well','Volume'])
                output = output.append(temp)
        return output.reset_index(drop=True)
            
class AssayPlate():
    def __init__(self):
        self.blocks = {}
        self.Allocations = {} # IDs map to blocks
        self.TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        
    def AddBlocks(self,block):
        count = len(self.blocks)+1
        self.blocks[count] = block

    def MapWells(self):
        # Gets the block from the number and replaces the
        # place holder destingation well IDs with actual well numbers
        TransferPlan = pd.DataFrame([],columns = ['SrcWell','Destination Well','Volume'])
        alphabet = string.ascii_uppercase
        for BlockNumber in self.blocks:
            # if odd
            if BlockNumber%2 == 1: 
                TestWells = [i+str(BlockNumber) for i in list(alphabet)[0:8]]
                BlankWells = [i+str(BlockNumber+1) for i in list(alphabet)[0:8]]
            else:
                TestWells = [i+str(BlockNumber-1) for i in list(alphabet)[8:16]]
                BlankWells = [i+str(BlockNumber) for i in list(alphabet)[8:16]]
                

            mappings = dict(zip(['X'+str(i) for i in range(1,9)]+['Y'+str(i) for i in range(1,9)],\
            TestWells + BlankWells))
            transfers = self.blocks[BlockNumber].Transfers
            transfers['Destination Well'] = transfers['Destination Well'].replace(mappings)
            self.TransferPlan = self.TransferPlan.append(transfers)
            self.TransferPlan.reset_index(inplace=True,drop=True)
    
    
   
def test_1():
    aracadonic = SourcePlateCompound('Arachadonic acid',['A'+str(i) for i in range(1,8)],ldv=False)
    DMSO = SourcePlateCompound('DMSO',['B'+str(i) for i in range(1,15)],ldv=False)
    assayplate = AssayPlate()
    for vol in [20,30,40]:
        for repeat in range(8):
            block = Block(aracadonic,DMSO,vol)
            assayplate.AddBlocks(block)
    assayplate.MapWells()
    transferMap = assayplate.TransferPlan
    transferMap = transferMap.loc[transferMap['Volume']>0]
    
    for i in transferMap['Destination Well'].unique():
        print(i, len(transferMap.loc[transferMap['Destination Well'] == i]))

def test_2():
    aracadonic = SourcePlateCompound('Arachadonic acid',['A'+str(i) for i in range(1,20)],ldv=False)
    DMSO = SourcePlateCompound('DMSO',['B'+str(i) for i in range(1,20)],ldv=False)

    assayplate1 = AssayPlate()
    assayplate2 = AssayPlate()
    assayplate3 = AssayPlate()
    assayplate4 = AssayPlate()
    
    for vol in [20,30,40]:
        for repeat in range(8):
            block1 = Block(aracadonic,DMSO,vol)
            assayplate1.AddBlocks(block1)
            block2 = Block(aracadonic,DMSO,vol)
            assayplate2.AddBlocks(block2)
            block3 = Block(aracadonic,DMSO,vol)
            assayplate3.AddBlocks(block3)
            block4 = Block(aracadonic,DMSO,vol)
            assayplate4.AddBlocks(block4)

    assayplate1.MapWells()
    assayplate2.MapWells()
    assayplate3.MapWells()
    assayplate4.MapWells()

    transferMap1 = assayplate1.TransferPlan
    transferMap2 = assayplate2.TransferPlan
    transferMap3 = assayplate3.TransferPlan
    transferMap4 = assayplate4.TransferPlan
    
    
if __name__ == '__main__':
    test_2()

import argparse
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters

def pick(smiles, n):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in mols]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    picker = SimDivFilters.MaxMinPicker()
    idx = picker.LazyPick(fn, len(smiles), n)
    return [smiles[i] for i in idx]

def lookup(smiles, df):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main(n):
    df = pd.read_csv(args.input)
    selection = lookup(pick(df['SMILES'], int(args.number)), df)
    selection.to_csv(args.output, 
            index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input')
    parser.add_argument('-o','--output')
    parser.add_argument('-n', '--number')
    args = parser.parse_args()
    main(args)
from string import ascii_uppercase
import pandas as pd
import echo

def main():
    wells = [f'{i}{j}' for i in ascii_uppercase[:16] for  j in range(1,25)]
    df = pd.read_csv('round1-cpds.csv')

    # dmso
    dmso = echo.Cpd(name = 'dmso', vol = 10_000)
    dmso_src = echo.Src(ldv = False)
    for i in wells[:24]:
        dmso_src.fill(i, dmso, 60)
    
    # compounds
    cpds = [echo.Cpd(name = i, vol = 20) for i in df['Item Name']]
    cpd_src = echo.Src()
    for i, j in zip(wells, cpds):
        cpd_src.fill(i, j)

    # blocks and dest fill
    dest = echo.Dest()
    # block1 = echo.Block(cpds[0], dmso, dest.wells) # slice?


if __name__ == '__main__':
    main()
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem import AllChem as Chem 
from rdkit.Chem import Draw

df = pd.read_csv('selection.csv', index_col=0).dropna()
mols = [Chem.MolFromSmiles(i) for i in df['SMILES']]
### pic
im=Draw.MolsToGridImage(mols, legends = df['Item Name'].tolist(), molsPerRow=10)
im.save('compound_selection_structures.png')
### hydrogens for property calcs
mols = [Chem.AddHs(i) for i in mols]


import numpy as np
import pandas as pd
import echo

# cpd needs to know which wells are its children 
# and send xfer command
# ? where to register xfer ?

def xfer_block(wells, a, b, vol = 1.5, k = 3):
    # wells: iterable of echo.Wells
    # a, b echo.Cpd
    # vol µl
    # k skew factor
    vols_a = (np.linspace(0,1,len(wells)) ** k) * vol
    vols_b = vol - vols_a
    for i,j,k in zip(wells, vols_a, vols_b):
        a.xfer(i,j)
        b.xfer(i,k)



def main():
    selected_compounds = pd.read_csv('selection.csv', index_col=0).dropna()

    cpds = [echo.Cpd(i, vol=100) for i in selected_compounds.index]
    dmso = echo.Cpd('dmso', vol=1500)

    src1 = echo.SrcPlate(name='src1')
    # next free well: J1 (n=211)
    src_dmso = echo.SrcPlate(name='src_dmso', ldv=False)

    src_plate_fill = []
    layout = []
    
    # used up to inc f24
    #In [13]: echo.hwells.index('F24')
    #Out[13]: 143
    for i,j,k in zip(cpds, src1[144::2], src1[145::2]):
        j.fill(i.sample(12))
        k.fill(i.sample(12))
        src_plate_fill.append([j.loc,i.name])
        src_plate_fill.append([k.loc,i.name])
        layout.append({i.name:[j.loc,k.loc]})

    for i in src_dmso[211:211+19]:
        i.fill(dmso.sample(60))

    for plate in ['bsa-blank', 'bsa-test','ctrl-blank','ctrl-test']:
        dest = echo.DestPlate(name=f'dest-{plate}')
        for i, j in zip(range(0,24*8,8), cpds):
            xfer_block(dest[i:i+8], j,dmso)

    src_xfers = pd.DataFrame(src1.xfer_record)
    dmso_xfers = pd.DataFrame(src_dmso.xfer_record)
    src_plate_fill = pd.DataFrame(src_plate_fill)
    src_plate_fill.to_csv('src_plate_fill.csv',index=False)
    src_xfers.to_csv('cpd-src-picklist.csv')
    dmso_xfers.to_csv('dmso-src-picklist.csv')

    ### layout suff
    layout = pd.DataFrame(layout).melt().dropna()
    layout.columns=['CatalogNumber','Src Wells']
    layout.index=layout['CatalogNumber']
    fda = pd.read_csv('fda-lib-layout.csv', index_col=0)
    fda.index = fda['CatalogNumber']
    fda = fda.reindex(layout['CatalogNumber'])
    layout = pd.concat([fda, layout['Src Wells']], axis=1)
    layout = layout[['Item Name', 'CatalogNumber', 'Rack Number',  'Plate Location',  'Src Wells']]  
    layout.to_csv('src_layout.csv', index=False)



if __name__ == '__main__':
    main()
import plates

def main():
    pass

if __name__ == '__main__':
    main()
import numpy as np
import pandas as pd
import echo

def make_block(wells, a, b, vol = 1.5, k = 3):
    # a nd b are compounds
    # vol is in µl
    # k is a skewing factor
    # return modified wells
    x = (np.arange(0,1, len(wells)) ** k) * vol
    for i,j in zip(wells, x):
         # change to xfer
         # required: parent plate/well, child plate/well
        i.fill(a, j)
        i.fill(b, vol - j)
    return wells


def main():
    df = pd.read_csv('selection.csv',index_col=0)
    dmso = echo.Cpd(name='dmso', vol = 1000)
    cpds = [echo.Cpd(name=i, vol=20) for i in df['CatalogNumber']]
    src_cpds = echo.Src()
    src_dmso = echo.Src(ldv=False)

    dest = echo.Dest()

    wells = make_block([dest.wells[i] for i in echo.hwells[:8]], cpds[0], dmso)
    print(wells)
    for i in wells:
        print(i.contents)

    # for i,j in zip(echo.hwells, cpds):
    #     src_cpds.fill(i, j, vol=12)

    # for i in echo.hwells[:5]:
    #     src_dmso.fill(i,dmso, 60)
    # 
    # for i in [1,2,3]:
    #     dest = echo.Dest()
    #     # make blocks


if __name__ == '__main__':
    main()
import pandas as pd
import cpds

def get_by_smiles(df, smiles):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main():
    df = pd.read_csv('fda-lib-layout.csv', index_col=0)
    smiles = df['SMILES']
    selection_smiles = cpds.pick(smiles, 24)

    selection = get_by_smiles(df, selection_smiles)
    selection.to_csv('selection.csv')

if __name__ == '__main__':
    main()
import numpy as np
import pandas as pd
import echo

# cpd needs to know which wells are its children 
# and send xfer command
# ? where to register xfer ?

def xfer_block(wells, a, b, vol = 1.5, k = 3):
    # wells: iterable of echo.Wells
    # a, b echo.Cpd
    # vol µl
    # k skew factor
    vols_a = (np.linspace(0,1,len(wells)) ** k) * vol
    vols_b = vol - vols_a
    for i,j,k in zip(wells, vols_a, vols_b):
        a.xfer(i,j)
        b.xfer(i,k)



def main():
    selected_compounds = pd.read_csv('selected_compounds.csv', index_col=0)

    cpds = [echo.Cpd(i, vol=100) for i in selected_compounds['CatalogNumber']]
    dmso = echo.Cpd('dmso', vol=1500)

    src1 = echo.SrcPlate(name='src1')
    # next free well: I1 (n=192)
    src_dmso = echo.SrcPlate(name='src_dmso', ldv=False)

    src_plate_fill = []
    
    for i,j,k in zip(cpds, src1[::2], src1[1::2]):
        j.fill(i.sample(12))
        k.fill(i.sample(12))
        src_plate_fill.append([j.loc,i.name])
        src_plate_fill.append([k.loc,i.name])

    for i in src_dmso[192:192+19]:
        i.fill(dmso.sample(60))

    for plate in ['bsa-blank', 'bsa-test','ctrl-blank','ctrl-test']:
        dest = echo.DestPlate(name=f'dest-{plate}')
        for i, j in zip(range(0,24*8,8), cpds):
            xfer_block(dest[i:i+8], j,dmso)
        dest.map.to_csv(f'{plate}-map.csv')

    src_xfers = src1.xfer_record
    dmso_xfers = src_dmso.xfer_record
    src_xfers.to_csv('cpd-src-picklist.csv')
    dmso_xfers.to_csv('dmso-src-picklist.csv')


if __name__ == '__main__':
    main()
import pandas as pd
import cpds

def main():
    df = pd.read_csv('fda-lib.csv', index_col = 0)
    selected_smiles = cpds.cpds.pick(df.SMILES, 25)

    selected_df = pd.concat([df.loc[df.SMILES == i,:] for i in selected_smiles])

    selected_df.to_csv('selected_compounds.csv',index=False)

if __name__ == '__main__':
    main()
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
import pandas as pd

def main():
    df = pd.read_excel('selleckchem-plate.xlsx',
            sheet_name = 'L1300-FDA-978cpds')
    df = df.loc[:,['Item Name','CatalogNumber', 'SMILES', 'Rack Number', 'Plate Location']]
    df.to_csv('layouts.csv')

if __name__ == '__main__':
    main()
import argparse
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters

def pick(smiles, n):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    fps = [Chem.RDKFingerprint(i) for i in mols]
    fn = lambda i, j : 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
    picker = SimDivFilters.MaxMinPicker()
    idx = picker.LazyPick(fn, len(smiles), n)
    return [smiles[i] for i in idx]

def lookup(smiles, df):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles])

def main(n):
    df = pd.read_csv(args.input)
    selection = lookup(pick(df['SMILES'], int(args.number)), df)
    selection.to_csv(args.output, 
            index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input')
    parser.add_argument('-o','--output')
    parser.add_argument('-n', '--number')
    args = parser.parse_args()
    main(args)
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen
from rdkit import DataStructs
from rdkit import SimDivFilters


def filter_HL(smiles):
    # filter smiles by herbicide likeness rules Hao
    mw_cutoff_min = 100
    mw_cutoff_max = 435
    logp_cutoff = 6
    hba_cutoff = 6
    hbd_cutoff = 2
    rotatable_cutoff =  9
    n_aromatic_bonds_cutoff = 17
    
    makemols = lambda smiles : Chem.AddHs(Chem.MolFromSmiles(smiles))
    n_aromatic_bonds = lambda m : sum([b.GetIsAromatic() for b in m.GetBonds()])
    mols = [makemols(i) for i in smiles]
    props = {s:{'mw':Chem.CalcExactMolWt(i), 
                'logp': Crippen.MolLogP(i),
                'hba': Chem.CalcNumLipinskiHBA(i),
                'hbd': Chem.CalcNumLipinskiHBD(i),
                'rot': Chem.CalcNumRotatableBonds(i),
                'aroB': n_aromatic_bonds(i)}
                for i,s in zip(mols, smiles)}
    
    prop_filter = lambda s : props[s]['mw'] <= mw_cutoff_max \
                        and props[s]['mw'] >= mw_cutoff_min \
                        and props[s]['logp'] <= logp_cutoff \
                        and props[s]['hba'] <= hba_cutoff\
                        and props[s]['hbd'] <= hbd_cutoff\
                        and props[s]['rot'] <= rotatable_cutoff\
                        and props[s]['aroB'] <= n_aromatic_bonds_cutoff
    return [i for i in smiles if prop_filter(i)]

def lookup(smiles, df):
    return pd.concat([df.loc[df['SMILES'] == i,:] for i in smiles]).drop('Unnamed: 0', axis = 1)

def main():
    df = pd.read_csv('layouts.csv')
    herbicideLike = filter_HL(df['SMILES'])
    herbicideLike = lookup(herbicideLike, df)
    # clean smiles
    herbicideLike['SMILES'] = [Chem.MolToSmiles(Chem.MolFromSmiles(i)) for i in herbicideLike['SMILES']]
    herbicideLike.to_csv('herbicide-like.csv', 
            index = False)

if __name__ == '__main__':
    main()
