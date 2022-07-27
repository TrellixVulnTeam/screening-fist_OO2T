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
from scipy.ndimage import convolve1d

import sxfst

'''
Data analysis script;
takes an input --root (directory) containing preprocessed traces
in csv(s) a la data_proc2.py.

outputs metrics and optionally plots too.
'''

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
    if isinstance(prot, float):
        x = df.loc[df['protein'].isna(), :]
    elif isinstance(prot, str):
        x = df.loc[df['protein'] == prot, :]
    x = x.loc[x['Cpd'].isna(), :]
    #well_row = sample['Well'].str.extract('([A-Z])')[0].unique()
    #assert len(well_row) == 1, f'{well_row}'
    #well_row = well_row[0]
    #test_run_no = sample['test_run_no'].unique()
    #test_run_no = test_run_no[0]
    #x = x.loc[x['Well'].str.contains(well_row), ]
    return x


def plotTraces(x,     # df
               ax,
               save_path=None,
               title=None,   # name
               concs=None,  #
               save=False,
               size=(12,8),
               ylim=None,#(-0.1,0.5),
               **kwargs,
               ):
    if concs is not None:
        for row_, conc_ in zip(x.index, concs):
            ax.plot(x.loc[row_,:],
                     c=plt.cm.inferno(conc_/max(concs)),
                     label=f'{round(conc_,2)} uM')
    else:
        for row_ in x.index:
            ax.plot(x.loc[row_,:],
                     label=f'{row_}',
                     **kwargs,
                     )
    ax.set_xlim(280,800)
    if ylim is not None:
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

def scale(x_, return_min_max=False):
    x = x_.copy()
    x -= min(x)
    x /= max(x)
    if return_min_max:
        return x, min(x_), max(x_)
    else:
        return x


def get_mm(x,y):
    x = np.nan_to_num(x, nan=1e-9) / 500
    y = np.nan_to_num(y, nan=1e-9)
    #xs, xmin, xmax = scale(x, return_min_max=True)
    #ys, ymin, ymax = scale(y, return_min_max=True)
    try:
        (km, vmax), covariance = curve_fit(mm, x, y,
                                           bounds=((1/32, max(y)*2), 
                                                   (4, max(y)*4)),
                                           p0=(1/8, max(y)*2),
                                           method='dogbox',
                                           #sigma=[1e-3]*len(x),
                                           )
        #(km_, vmax_), covariance = curve_fit(mm, xs, ys,
        #                                   bounds=((min(xs), min(ys)), 
        #                                           (max(xs), max(ys))),
        #                                   )
        #km = (km_ * ymax) + ymin
        #vmax = (vmax_ * xmax) + xmin
        km *= 500
        x  *= 500
    except RuntimeError:
        km, vmax = np.inf, np.inf

    yh = mm(x, km=km, vmax=vmax)
    rsq = r_squared(y, yh)
    return {'km':km, 'vmax':vmax, 'rsq':rsq}

def get_extra_metrics(test_traces, ctrl_traces):
    def has420peak(traces):
        pass
    return {}


def plot_mm(ax,
            x,
            y,
            km,
            vmax,
            rsq,
            title,
           ):
    ax.scatter(x,y)
    xx = np.linspace(min(x),max(x), 64)
    yy = mm(xx, vmax, km)
    ax.plot(xx, yy)
    ax.set_xlabel('Concenctration uM')
    ax.set_ylabel('Response')
    ax.set_title(title)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlim(-8, 512)
    ax.text(x=350,
            y=0.9,
            s=f'kd: {round(km, 5)}\nvmax: {round(vmax, 5)}\nrsq: {round(rsq,3)}',
            fontsize=14,
            )
    ax.hlines(y=vmax, 
              xmin=-8,
              xmax=max(x),
              linestyle='--',
              lw=1,
              color='gray',
              )
    ax.hlines(y=vmax/2, 
              xmin=-8,
              xmax=km,
              linestyle='--',
              lw=1,
              color='gray',
              )
    ax.vlines(x=km, 
              ymin=-0.1,
              ymax=mm(km, vmax, km),
              linestyle='--',
              lw=1,
              color='gray',
              )


def trace_similarity(a, b):
    similarity = sorted(range(len(protein_blanks_traces)),
                        key=lambda idx : abs(protein_blanks_traces.iloc[idx, 400].mean() - \
                                        test_traces.iloc[:,400].mean()))
    pass


def main(args):
    plt.style.use('default')
    root = args.data
    img_root = args.img
    if args.out != '':
        if not os.path.exists(args.out):
            os.mkdir(args.out)
    csvs = [os.path.join(root, i) for i in os.listdir(root)]
    df = pd.concat([pd.read_csv(i, low_memory=False) for i in csvs]).reset_index(drop=True)
    if img_root is not None:
        img_paths = [i for i in sxfst.find(img_root) if 'png' in i]  ###

    header_done = False # True if csv header already written
    sigma = 2
    for i in df['protein'].dropna().unique():
        for j in tqdm(df['Cpd'].dropna().unique(), disable=args.stdout):
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
                                    key=lambda idx : abs(protein_blanks_traces.iloc[idx, 400].mean() - \
                                                    test_traces.iloc[:,400].mean()))
                protein_blanks_trace = sxfst.data.smooth(protein_blanks_traces.iloc[[similarity[0]],:], 
                                                         sigma=sigma)
                
                test_traces = pd.concat([protein_blanks_trace,
                                          test_traces],
                                        axis=0)
                
                control_blanks = get_blank_wells(df, ctrl)
                control_blanks_traces = get_traces(control_blanks)
                similarity = sorted(range(len(control_blanks_traces)),
                                    key=lambda idx : abs(control_blanks_traces.iloc[idx, 300].mean() - \
                                                    test_traces.iloc[:,300].mean()))
                control_blanks_trace = sxfst.data.smooth(control_blanks_traces.iloc[similarity[0],:],
                                                         sigma=sigma)#Series
                control_blanks_trace = control_blanks_trace.sub(control_blanks_trace.iloc[-1].values, axis=1)
                ctrl_traces = get_traces(ctrl)
                ctrl_traces = pd.concat([pd.DataFrame(control_blanks_trace).T,
                                          get_traces(ctrl)],
                                          axis=0)
                ctrl_traces_norm = sxfst.data.norm_traces(ctrl_traces)
                ctrl_traces_norm = ctrl_traces_norm.sub(ctrl_traces_norm.iloc[0,:].values)
                assert sum(ctrl_traces_norm.loc[:,800]) == 0 , f'{ctrl_traces_norm}'
                ctrl_traces_norm_sub = ctrl_traces_norm.sub(control_blanks_trace[0], axis=1)
                #.sub(control_blanks_trace.values, axis=1)
                #ctrl_traces_norm_ = sxfst.data.norm_traces(ctrl_traces)
                #ctrl_traces_norm = ctrl_traces_norm_.sub(control_blanks_trace, 
                #                                         axis=1)
                ctrl_traces_smooth = sxfst.data.smooth(ctrl_traces_norm, sigma=sigma)
                #print(ctrl_traces_norm)
                #print(ctrl_traces_smooth)

                #vols = test['actual_vol'].to_list()
                vols = [0] + test['actual_vol'].to_list()
                concs = np.array([sxfst.data.c2(v1=i,      # vol
                                                c1=10_000, # stock conc - uM
                                                v2=38_000 + i, # total vol nm
                                                ) for i in vols])

                test_traces_norm = sxfst.data.norm_traces(test_traces)
                #test_traces_norm_ = sxfst.data.norm_traces(test_traces)
                #test_traces_norm = test_traces_norm_.sub(control_blanks_trace, 
                #                                         axis=1)
                test_traces_smooth = sxfst.data.smooth(test_traces_norm, 
                                                        sigma=sigma, 
                                                        axis=1).sub(control_blanks_trace.iloc[:,0], 
                                                                    axis=1)
                def gradient(df):
                    x = convolve1d(test_traces_smooth, [-1,0,1])
                    return pd.DataFrame(x, columns=df.columns, index=df.index)

                grad = gradient(test_traces_smooth)
                diff = grad - grad.iloc[0,:]
                #diff = test_traces_smooth - test_traces_smooth.iloc[0,:]
                
                response = sxfst.data.response(grad.sub(grad.iloc[0,:].values), a=410, b=439)
                
                mm_fit = get_mm(concs, response.values) # dict
                
                extra_metrics = get_extra_metrics(test_traces_smooth, ctrl_traces) # dict

                output_data = {'cpd': j,
                               'protein' : i,
                               **mm_fit,
                               **extra_metrics,
                               }
                odf = pd.DataFrame({0:output_data}).T
                odf = odf.loc[:,['cpd','protein','km','vmax','rsq']]

                if args.stdout:
                    if not header_done:
                        odf.to_csv(sys.stdout, index=False)
                        header_done = True
                    else:
                        odf.to_csv(sys.stdout, index=False, header=False)
                elif args.file_out:
                    if not header_done:
                        odf.to_csv(os.path.join(args.out, args.file_out), index=False)
                        header_done = True
                    else:
                        odf.to_csv(os.path.join(args.out, args.file_out), mode='a',index=False, header=False)
                if args.plot:
                    fig, ax = plt.subplots(3,2, figsize=(16, 16))
                    plotTraces(test_traces_smooth,
                               concs=concs,
                               ylim=(-0.05, 0.25),
                               size=(8,3),
                               ax=ax[0,0],
                               title=f'{i} : {j} - Test Traces',
                               )
                    ax[0,0].hlines(y=0, 
                                   xmin=0, 
                                   xmax=800, 
                                   linestyle='--',
                                   lw=1,
                                   color='gray',
                                   )
                    ax[0,0].vlines([390, 420], 
                                   [0, 0], 
                                   [test_traces_smooth[390].max(), test_traces_smooth[420].max()], 
                                   linestyle='--',
                                   lw=1,
                                   color='gray',
                                   )
                    plotTraces(ctrl_traces_smooth,
                               concs=concs,
                               ylim=(-0.05,0.25),
                               size=(8,3),
                               ax=ax[0,1],
                               title=f'{i} : {j} - Control Traces',
                               )
                    ax[0,1].hlines(y=0, 
                                   xmin=0, 
                                   xmax=800, 
                                   linestyle='--',
                                   lw=1,
                                   color='gray',
                                   )
                    plotTraces(grad,
                               concs=concs,
                               ylim=(-0.02, 0.03),
                               size=(8,3),
                               ax=ax[1,0],
                               title=f'{i} : {j} - Trace Gradients',
                               )
                    ax[1,0].hlines(y=0, 
                                   xmin=0, 
                                   xmax=max(concs), 
                                   linestyle='--',
                                   lw=1,
                                   color='gray',
                                   )
                    ax[1,0].vlines([410, 430], 
                                   [0, 0], 
                                   [grad[410].min(), grad[430].max()], 
                                   linestyle='--',
                                   lw=1,
                                   color='gray',
                                   )

                    plot_mm(ax[1,1],
                            x=concs,
                            y=response,
                            km=mm_fit['km'],
                            vmax=mm_fit['vmax'],
                            rsq=mm_fit['rsq'],
                            title=f"{i} : {j} - Michaelis Menten Fit",
                           )
                    ax[1,1].set_ylim=(-1e-3, 1e-3)

                    if img_root is not None:
                        cpd_img = Image.open(next(filter(lambda s : j in s, img_paths)))
                        ax[2,0].imshow(cpd_img)
                        ax[2,0].axis('off')
                    ax[2,1].axis('off')

                    plt.tight_layout()

                    save_path = os.path.join(args.out, 
                            f'{i}:{j}.png'.replace(' ','-').replace('/', '_'))
                    plt.savefig(save_path)
                    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', help='root data dir')
    parser.add_argument('-i', '--img', help='root img dir')
    parser.add_argument('-p', '--plot', help='plot (bool)', action='store_true')
    parser.add_argument('-o', '--out', default='', help='output dir name')
    parser.add_argument('-s', '--stdout', help='write to stdout', action='store_true')
    parser.add_argument('-f', '--file_out', help='out csv name, default=out.csv', default='out.csv')
    args = parser.parse_args()
    main(args)
