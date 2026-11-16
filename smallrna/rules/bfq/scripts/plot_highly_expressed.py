#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import pandas as pd
import numpy as np
import scanpy as sc

from scipy.sparse import issparse

def highly_expressed_yaml(adata, n_top=10, biotype=None):
    """The n_top miRNAs with the highest mean counts per million miRNA reads (CPM), as CPM and as read counts.
    """
    adata = adata.copy()
    norm_dict = sc.pp.normalize_total(adata, target_sum=1e6, inplace=False)
    adata.layers['CPM'] = norm_dict['X']
    adata.obs['norm_factor'] =  norm_dict['norm_factor']
    
    if issparse(norm_dict['X']):
        mean_cpm = norm_dict['X'].mean(axis=0).A1
        top_idx = np.argsort(mean_cpm)[::-1][:n_top]
    else:
        mean_cpm = norm_dict['X'].mean(axis=0)
        top_idx = np.argsort(mean_cpm)[::-1][:n_top]
    top_genes = adata.var_names[top_idx]
    top_adata = adata[:,top_genes]
    
    keys = {}
    default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c',
                      '#8085e9','#f15c80', '#e4d354', '#2b908f',
                      '#f45b5b', '#91e8e1']
    for i, name in enumerate(top_genes):
        color = default_colors[i]
        keys[name] = {'color': color, 'name': name}
    
    # Config for the plot
    data_labels= [{'name': 'CPM', 'ylab': 'CPM'}, {'name': 'counts', 'ylab': '# reads'}]
    data = [top_adata.to_df(layer='CPM').to_dict(orient='index'), top_adata.to_df().to_dict(orient='index')]
    pconfig = {
        'id': 'gene_high',
        'cpswitch': False,
        'title': 'Highly expressed miRNAs',
        'data_labels': data_labels
    }
    
    section = {}
    section['id'] = 'gene_high'
    section['section_name'] = 'Highly expressed miRNAs'
    section['description'] = ('The {} miRNAs with the highest mean share of miRNA reads across samples, as counts per million '
                              'miRNA reads in each sample (CPM) or as read counts.'.format(n_top))
    section['plot_type'] = 'bargraph'
    section['pconfig'] = pconfig
    section['categories'] = keys
    section['data'] = data

    return section


def argparser():
    parser = argparse.ArgumentParser(description="Highly expressed miRNAs figure for QC report")
    parser.add_argument("adata", help="scanpy formatted file (.h5ad0) ")
    parser.add_argument("-o ", "--output", default="top10_mqc.pyaml", help="Output filename. Will default to top10_mqc.yaml")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = argparser()
    adata = sc.read(args.adata)
    section =  highly_expressed_yaml(adata)
    with open(args.output, 'w') as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
