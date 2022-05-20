#!/usr/bin/env python
import sys
import os
import argparse
import json
from tqdm import tqdm

from PIL import Image
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import sxfst

def get_experiment(df, 
                   protein, 
                   compound):
    test = df.loc[df['Cpd'] == compound,:].loc[df['protein'] == protein,:]
    ctrl = df.loc[df['Cpd'] == compound,:].loc[df['protein'].isna(),:]
    return test, ctrl

def get_traces(df):
    traces = df.loc[:,'300':]
    traces.columns = traces.columns.astype(int)
    traces.index = df['Well']
    return traces

def get_blank_wells(df, sample):
    prot = sample['protein'].unique()
    assert len(prot) == 1, f'{prot}'
    prot = prot[0]
    
    well_row = sample['Well'].str.extract('([A-Z])')[0].unique()
    assert len(well_row) == 1, f'{well_row}'
    well_row = well_row[0]
    
    test_run_no = sample['test_run_no'].unique()
    test_run_no = test_run_no[0]
    x = df.loc[df['test_run_no'] == test_run_no, :]
    x = x.loc[x['Cpd'].isna(), :]
    x = x.loc[x['Well'].str.contains(well_row), ]
    return x


def plotTraces(x,     # df
               ax,
               save_path=None,
               title=None,   # name
               concs=None,  #
               save=False,
               size=(12,8),
               ylim=(-0.1,0.5),
               **kwargs,
               ):
    if concs is not None:
        for row_, conc_ in zip(x.index, concs):
            ax.plot(x.loc[row_,:],
                     c=plt.cm.cool(conc_/max(concs)),
                     label=f'{round(conc_,2)} uM')
    else:
        for row_ in x.index:
            ax.plot(x.loc[row_,:],
                     label=f'{row_}',
                     **kwargs,
                     )
    ax.set_xlim(280,800)
    ax.set_ylim(*ylim)
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Absorbance')
    ax.legend(loc='center right')
    return ax

def mm(x, vmax, km):
    return ((x * vmax) / (km + x))


def r_squared(yi,yj):
    residuals = yi - yj
    sum_sq_residual = sum(residuals ** 2)
    sum_sq_total = sum((yi - yi.mean()) ** 2) # check this!!!
    return 1 - (sum_sq_residual / sum_sq_total)

def get_mm(x,y):
    #y = y.replace(np.inf, 0) # error handling - pandas
    try:
        (km, vmax), covariance = curve_fit(mm, x, y,
                bounds=((0, 0),(1e4,0.5)),
                p0 = (1e4, 0))
    except RuntimeError:
        km, vmax = np.inf, np.inf

    yh = mm(x, km, vmax)
    rsq = r_squared(y, yh)
    return {'km':round(km,2), 'vmax':round(vmax,2), 'rsq':round(rsq,2)}


def plot_mm(ax,
            x,
            y,
            km,
            vmax,
            title,
           ):
    ax.scatter(x,y)
    xx = np.linspace(min(x),max(x), 64)
    yy = mm(xx, vmax, km)
    ax.plot(xx, yy)
    ax.set_xlabel('Concenctration uM')
    ax.set_ylabel('Response')
    ax.set_title(title)




def main(args):
    root = args.data
    img_root = args.img
    odir = 'out'
    if not os.path.exists(odir):
        os.mkdir(odir)
    csvs = [os.path.join(root, i) for i in os.listdir(root)]
    df = pd.concat([pd.read_csv(i, low_memory=False) for i in csvs]).reset_index(drop=True)
    img_paths = [i for i in sxfst.find(img_root) if 'png' in i]  ###

    header_done = False # True if csv header already written
    for i in df['protein'].dropna().unique():
        if args.stdout:
            tqdm = lambda x : x
        for j in tqdm(df['Cpd'].dropna().unique()):
            test, ctrl = get_experiment(df, i, j)
            test_run_no = test['test_run_no'].unique()
            assert len(test_run_no) == 1, f'{test_run_no, i, j }'
            test_run_no = test_run_no[0]
            ctrl = ctrl.loc[ctrl['test_run_no'] == test_run_no, :]
            
            if len(test) > 0:
                test_traces = get_traces(test)
                protein_blanks = get_blank_wells(df, test)
                protein_blanks_traces = get_traces(protein_blanks)
                # get most similar at A400 - come back to this
                similarity = sorted(range(len(protein_blanks_traces)),
                                    key=lambda idx : abs(protein_blanks_traces.iloc[idx, 400].mean() - test_traces.iloc[:,400].mean()))
                protein_blanks_trace = protein_blanks_traces.iloc[[similarity[0]],:]
                
                test_traces = pd.concat([protein_blanks_trace,
                                         test_traces],
                                       axis=0)
                
                control_blanks = get_blank_wells(df, ctrl)
                control_blanks_traces = get_traces(control_blanks)
                control_blanks_trace = control_blanks_traces.loc[[control_blanks_traces.index[0]],:]
                ctrl_traces = pd.concat([control_blanks_trace,
                                         get_traces(ctrl)],
                                         axis=0)
                vols = [0] + test['actual_vol'].to_list()
                concs = np.array([sxfst.data.c2(v1=i,      # vol
                                                c1=10_000, # stock conc - uM
                                                v2=38 + i, # total vol
                                                )
                                     for i in vols])

                test_traces_norm = sxfst.data.norm_traces(test_traces)
                test_traces_smooth = sxfst.data.smooth(test_traces_norm)
                diff = test_traces_smooth - test_traces_smooth.iloc[0,:]
                response = sxfst.data.response(diff)
                mm_fit = get_mm(np.array(concs), response.values)
                
                output_data = {'cpd': j,
                               'protein' : i,
                               **mm_fit}
                odf = pd.DataFrame({0:output_data}).T

                if args.stdout:
                    odf.to_csv(sys.stdout)
                elif args.file_out:
                    if not header_done:
                        odf.to_csv(args.file_out, index=False)
                        header_done = True
                    else:
                        odf.to_csv(args.file_out, mode='a',index=False, header=False)
                if args.plot:
                    fig, ax = plt.subplots(2,2, figsize=(16,8))
                    plotTraces(test_traces_smooth,
                               concs=concs,
                               ylim=(-0.05, 0.3),
                               size=(8,3),
                               ax=ax[0,0],
                               title=f'{i} : {j} - Normalized Traces',
                               )
                    plotTraces(diff,
                               concs=concs,
                               ylim=(-0.1,0.1),
                               size=(8,3),
                               ax=ax[0,1],
                               title=f'{i} : {j} - Difference Traces',
                               )
                    
                    plot_mm(ax[1,0],
                            x=concs,
                            y=response,
                            km=mm_fit['km'],
                            vmax=mm_fit['vmax'],
                            title=f'{i} : {j} - Michaelis Menten - {mm_fit}',
                           )
                    cpd_img = Image.open(next(filter(lambda s : j in s,
                                                     img_paths)))
                    ax[1,1].imshow(cpd_img)
                    ax[1,1].axis('off')
                    plt.tight_layout()

                    save_path = os.path.join(odir, f'{i}:{j}.png'.replace(' ','-')
                    plt.savefig(save_path)
                    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', help='root data dir')
    parser.add_argument('-i', '--img', help='root img dir')
    parser.add_argument('-p', '--plot', help='plot (bool)', action='store_true')
    parser.add_argument('-o', '--out', help='output dir name')
    parser.add_argument('-s', '--stdout', help='write to stdout', action='store_true')
    parser.add_argument('-f', '--file_out', help='out csv name, default=out.csv', default='out.csv')
    args = parser.parse_args()
    main(args)
