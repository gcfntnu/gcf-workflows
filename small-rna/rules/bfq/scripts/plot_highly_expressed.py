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

def higly_expressed_yaml(adata, n_top=10, biotype=None):
    """plot the top 10 highest expressed genes
    """
    adata = adata.copy()
    norm_dict = sc.pp.normalize_total(adata, target_sum=100, inplace=False)
    adata.layers['CPM'] = norm_dict['X']
    adata.obs['norm_factor'] =  norm_dict['norm_factor']
    
    if issparse(norm_dict['X']):
        mean_percent = norm_dict['X'].mean(axis=0).A1
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
    else:
        mean_percent = norm_dict['X'].mean(axis=0)
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
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
        'title': 'Higly Expressed miRNA',
        'data_labels': data_labels
    }
    
    section = {}
    section['id'] = 'gene_high'
    section['section_name'] = 'Higly Expressed miRNA'
    section['description'] = 'Summary of the 10 most highly expressed miRNA.'
    section['plot_type'] = 'bargraph'
    section['pconfig'] = pconfig
    section['categories'] = keys
    section['data'] = data

    return section


def argparser():
    parser = argparse.ArgumentParser(description="Higly expressed genes figure for QC report")
    parser.add_argument("adata", help="scanpy formatted file (.h5ad0) ")
    parser.add_argument("-o ", "--output", default="top10_mqc.pyaml", help="Output filename. Will default to top10_mqc.yaml")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = argparser()
    adata = sc.read(args.adata)
    section =  higly_expressed_yaml(adata)
    with open(args.output, 'w') as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
