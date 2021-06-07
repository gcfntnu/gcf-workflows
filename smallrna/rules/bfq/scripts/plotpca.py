#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings('ignore', message='numpy.dtype size changed')

import pandas as pd
import numpy as np
import scanpy as sc
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')
sc.settings.autoshow = False
sc.settings.figdir = '.'
sc.plot_suffix = ''


import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

def argparser():
    parser = argparse.ArgumentParser(description='Dimred plot from anndata')
    parser.add_argument('-i', '--input',
                        help='Input filename. Anndata file (.h5ad)')
    parser.add_argument('-o', '--output', default='umap_mqc.png',
                        help='Output filename. Will default to umap_mqc.png, Optional')
    parser.add_argument('--recipe', default='mirna', choices = ['smallrna', 'rna', 'microbiome', 'deicode', 'vsd'],
                        help='Preprocess recipe. ')
    parser.add_argument('--method', default='auto', choices = ['auto', 'umap', 'pca'],
                        help='Dimension reduction method.')    
    args = parser.parse_args()
    return args


def multiqc_png(adata, args):
    if use_umap:
        if args.output.endswith('.png'):
            sc.pl.umap(adata, color=color, save=args.output)
        else:
            tab_pc = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])  
            tab_pc.index.name = 'Sample_ID'
            tab_pc.to_csv(args.output, sep='\t')
    else:
        if args.output.endswith('.png'):
            sc.pl.pca(adata, color=color, save=args.output)
        else:
            tab_pc = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names, columns=['PC1', 'PC2'])  
            tab_pc.index.name = 'Sample_ID'
            tab_pc.to_csv(args.output, sep='\t')
    

def multiqc_section(args):
    section =  {}
    section['parent_id'] = 'exploratory'
    section['parent_name'] = 'Exploratory analysis'
    section['parent_description'] = 'Exploratory analysis of feature table. A feature table is a count matrix with feature counts for each sample.'
    
    if args.method == 'pca':
        section['id'] = 'pca'
        section['section_name'] = 'PCA'
        section['description'] = 'Principal Component Analysis of tabulated counts.'
    elif args.method == 'umap':
        section['id'] = 'umap'
        section['section_name'] = 'UMAP'
        section['description'] = 'Uniform Manifold Approximation of tabulated counts.'  
    section['plot_type'] = 'scatter'
    return section

def multiqc_pconfig(args):
    pconfig = {}
    pconfig['id'] = 'pca_scatter'
    pconfig['title'] = 'PCA'
    pconfig['xlab'] = 'PC1'
    pconfig['ylab'] = 'PC2'
    pconfig['tt_label'] =  'PC1 {point.x:.2f}: PC2 {point.y:.2f}'
    return pconfig

def multiqc_yaml(adata, args):
    if 'X_umap' in adata.obsm:
        x = adata.obsm['X_umap'].copy()
    else:
        x = adata.obsm['X_pca'].copy()
    max_comp = x.shape[1]
    T = pd.DataFrame(x[:,:2], index=adata.obs_names, columns=['x', 'y'])  
    T.index.name = 'Sample_ID'
    S = adata.obs.copy()
    df1 = pd.concat([T, S], axis="columns")
    if max_comp > 3:
        T2 = pd.DataFrame(x[:,2:4], index=adata.obs_names, columns=['x', 'y'])  
        T2.index.name = 'Sample_ID'
        df2 = pd.concat([T2, S], axis="columns")
    
    section =  multiqc_section(args)
    pconfig = multiqc_pconfig(args)
    

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
        data1.update(T.to_dict(orient='index'))
    data.append(data1)
    
    if max_comp > 3:
        df2 = pd.concat([T2, S], axis="columns")
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
            data2.update(T2.to_dict(orient='index'))    
        data.append(data2)
        
    section['pconfig'] = pconfig
    section['data'] = data
    return section

def rpca(adata):
    from deicode.matrix_completion import MatrixCompletion
    from deicode.preprocessing import rclr
    min_samples = max(3, np.floor(n_samples * 0.1))
    sc.pp.filter_genes(adata, min_cells=min_samples) 
    X = rclr(adata.raw.X)
    opt = MatrixCompletion(n_components=n_comps, max_iterations=10).fit(X)
    n_components = opt.s.shape[0]
    X = opt.sample_weights @ opt.s @ opt.feature_weights.T
    X = X - X.mean(axis=0)
    X = X - X.mean(axis=1).reshape(-1, 1)
    adata.obsm['X_deicode'] = sc.tl.pca(X, svd_solver='arpack', n_comps=n_comps)
    return adata

if __name__ == "__main__":
    args = argparser()
    
    adata = sc.read(args.input)

    #adata.obs['Sample_Group'] = adata.obs['Tissue'].copy()
    
    n_samples, n_genes = adata.shape
    if args.method == 'auto':
        if n_samples >= 20:
            args.method = 'umap'
        else:
            args.method = 'pca'

    color = ['Sample_Group'] if 'Sample_Group' in adata.obs.columns else None
    adata.raw = adata
    n_comps = min(20, min(adata.shape)-1)
    
    if args.recipe == 'mirna':
        min_samples = max(3, np.floor(n_samples * 0.1))
        sc.pp.filter_genes(adata, min_cells=min_samples)
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.1, max_mean=20, min_disp=0.2)
        adata = adata[:, adata.var.highly_variable]
        
    elif args.recipe == 'deicode':
        min_samples = max(3, np.floor(n_samples * 0.1))
        sc.pp.filter_genes(adata, min_cells=min_samples) 
        adata = rpca(adata)
        
    elif args.recipe == 'mrna':
        sc.pp.filter_genes(adata, min_cells=min_samples)
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.8, max_mean=20, min_disp=0.2)
        adata = adata[:, adata.var.highly_variable]

    elif args.recipe == 'vsd':
        pass
        
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)

    if args.method == 'umap':
        sc.pp.neighbors(adata, n_neighbors=min(7, np.ceil(n_samples/3.0)), n_pcs=min(15, n_comps))
        sc.tl.leiden(adata, resolution=0.5)
        sc.tl.paga(adata)
        sc.pl.paga(adata, plot=False)
        sc.tl.umap(adata, init_pos='paga')

    if args.method == 'umap':
        if args.output.endswith('_mqc.png'):
            fig = sc.pl.umap(adata, color=color, show=False, return_fig=True)
            fig.savefig(args.output, dpi=300, transparent=True)
        elif args.output.endswith('_mqc.yaml'):
            section = multiqc_yaml(adata, args)
            with open(args.output, 'w') as fh:
                yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
        else:
            tab_pc = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])  
            tab_pc.index.name = 'Sample_ID'
            tab_pc.to_csv(args.output, sep='\t')
    else:
        if args.output.endswith('_mqc.png'):
            fig = sc.pl.pca(adata, color=color, return_fig=True, show=False)
            fig.savefig(args,output, dpi=300, transparent=True)
        elif args.output.endswith('_mqc.yaml'):
            section = multiqc_yaml(adata, args)
            with open(args.output, 'w') as fh:
                yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
        else:
            tab_pc = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names, columns=['PC1', 'PC2'])  
            tab_pc.index.name = 'Sample_ID'
            tab_pc.to_csv(args.output, sep='\t')

