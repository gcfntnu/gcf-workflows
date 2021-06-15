#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

import scanpy as sc


def var_filter(E, f=0.5):
    E = E.loc[E.sum(1) > 10, :]
    o = E.var(axis=1).argsort()[::-1]
    n = int(np.floor(E.shape[0] * f))
    E = E.iloc[o[:n], :]
    return E


def pca(E, max_comp=20):
    x = E.values
    max_comp = min(max_comp, min(x.shape))
    x = x - x.mean(0)
    u, s, vt = np.linalg.svd(x, 0)
    T = u * s
    P = vt.T
    T = pd.DataFrame(
        T[:, :max_comp],
        index=E.index,
        columns=["PC_{}".format(i + 1) for i in range(max_comp)],
    )
    P = pd.DataFrame(
        P[:, :max_comp],
        index=E.columns,
        columns=["PC_{}".format(i + 1) for i in range(max_comp)],
    )
    return T, P


def scanpy_data(E, S=None, F=None):
    sc.

def pairs_scores(T, S=None, n_comp=4):
    n_comp = min(T.shape[1], n_comp)
    T = T.iloc[:, :n_comp]
    df = pd.concat([T, S], axis="columns")
    sns.set(style="ticks")
    hue = "Sample_Group" if "Sample_Group" in df.columns else None
    ax = sns.pairplot(df, vars=T.columns, hue=hue)
    return ax

def multiqc_fig(T, S=None, n_comp=4):
    p = pairs_scores(T, S, n_comp=n_comp)
    plt.subplots_adjust(top=0.9)
    p.fig.suptitle("PCA on variance transformed gene expression data")
    p.savefig(args.output)
    
def multiqc_yaml(T, S=None, n_comp=2):
    max_comp = T.shape[1]
    
    T1 = T.iloc[:, :2]
    T1.columns = ['x', 'y']
    df1 = pd.concat([T1, S], axis="columns")   

    if max_comp > 3:
        T2 = T.iloc[:, 2:4]
        T2.columns = ['x', 'y']
        df2 = pd.concat([T2, S], axis="columns")   
        
    section =  {}
    section['id'] = 'pca'
    section['section_name'] = 'PCA'
    section['description'] = 'Principal Component Analysis of variance transformed gene counts.'
    section['plot_type'] = 'scatter'
    
    pconfig = {}
    pconfig['id'] = 'pca_scatter'
    pconfig['title'] = 'PCA'
    pconfig['xlab'] = 'PC1'
    pconfig['ylab'] = 'PC2'
    pconfig['tt_label'] =  'PC1 {point.x:.2f}: PC2 {point.y:.2f}'
    #pconfig['xmax'] = float(T['x'].max() + 0.1*T['x'].max())
    #pconfig['xmin'] = float(T['x'].min() + 0.1*T['x'].min())
    #pconfig['ymax'] = float(T['y'].max() + 0.1*T['y'].max())
    #pconfig['ymin'] = float(T['y'].min() + 0.1*T['y'].min())

    data_labels = [{'name': 'PC1 vs PC2', 'xlab': 'PC1', 'ylab': 'PC2'}]
    if max_comp > 3:
        data_labels.append({'name': 'PC3 vs PC4', 'xlab': 'PC3', 'ylab': 'PC4'})
    pconfig['data_labels'] = data_labels

    data = []
    data1 = {}
    default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9',
                      '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']
    if 'Sample_Group' in df1.columns:
        groups = set(df1['Sample_Group'])
        g_df = df1.groupby('Sample_Group')
        for i, g in enumerate(groups):
            sub = g_df.get_group(g)[['x', 'y']]
            sub['color'] = default_colors[i]
            sub['name'] = g
            this_series = sub.to_dict(orient='index')
            data1.update(this_series)
    else:
        data1.update(T1.to_dict(orient='index'))
    data.append(data1)
    
    if max_comp > 3:
        data2 = {}
        if 'Sample_Group' in df2.columns:
            groups = set(df2['Sample_Group'])
            g_df = df2.groupby('Sample_Group')
            for i, g in enumerate(groups):
                sub = g_df.get_group(g)[['x', 'y']]
                sub['color'] = default_colors[i]
                sub['name'] = g
                this_series = sub.to_dict(orient='index')
                data2.update(this_series)
        else:
            data2.update(T1.to_dict(orient='index'))    
        data.append(data2)
        
    section['pconfig'] = pconfig
    section['data'] = data
    
    
    with open(args.output, 'w') as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)

def argparser():
    parser = argparse.ArgumentParser(description="PCA figure for QC report")
    parser.add_argument("exprs")
    parser.add_argument(
        "--sample-info",
        help="Optional sample sheet. Will subset expr table if needed",
        dest="samples",
    )
    parser.add_argument(
        "--feature-info",
        help="Optional sample sheet. Will subset expr table if needed",
        dest="features",
    )
    parser.add_argument(
        "-o ",
        "--output",
        default="pca_mqc.png",
        help="Output filename. Will default to pca_mqc.png, Optional [*_mqc.yaml]",
    )
    parser.add_argument("--include-comp-test", action="store_true", dest="test_comp")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argparser()

    E = pd.read_csv(args.exprs, sep="\t", index_col=0)
    E.columns = E.columns.astype(str)
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not E.columns.isin(S.index).all():
            raise ValueError("missing samples in sample info!")
        S = S.loc[E.columns, :]
    if args.features is not None:
        F = pd.read_csv(args.features, sep="\t", index_col=0)
        if not E.index.isin(F.index).all():
            warnings.warn("missing annotations in feature info!")
        F = F.loc[E.index, :]

    F = var_filter(E)
    T, P = pca(F.T)

    if args.output.endswith('_mqc.png'):
        multiqc_fig(T, S)
    elif args.output.endswith('_mqc.yaml'):
        multiqc_yaml(T, S)
    else:
        raise ValueError('Output file need to end with _mqc.png or _mqc.yaml')
