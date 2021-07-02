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
