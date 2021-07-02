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
