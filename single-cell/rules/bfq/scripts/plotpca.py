#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import scanpy as sc
import pandas as pd

import matplotlib
matplotlib.use('Agg')

def multiqc_yaml(T, S=None, n_comp=2, color_by="Sample_Group"):
    max_comp = T.shape[1]
    T1 = T.iloc[:, :2]
    T1.columns = ['x', 'y']
    df1 = pd.concat([T1, S], axis="columns")   

    if max_comp > 3:
        T2 = T.iloc[:, 2:4]
        T2.columns = ['x', 'y']
        df2 = pd.concat([T2, S], axis="columns")   
        
    section =  {}
    section['id'] = 'umap'
    section['section_name'] = 'umap'
    section['description'] = 'UMAP Analysis of variance transformed gene counts.'
    section['plot_type'] = 'scatter'
    
    pconfig = {}
    pconfig['id'] = 'umap_scatter'
    pconfig['title'] = 'UMAP'
    pconfig['xlab'] = 'UMAP1'
    pconfig['ylab'] = 'UMAP2'
    pconfig['tt_label'] =  'UMAP1 {point.x:.2f}: UMAP2 {point.y:.2f}'
    #pconfig['xmax'] = float(T['x'].max() + 0.1*T['x'].max())
    #pconfig['xmin'] = float(T['x'].min() + 0.1*T['x'].min())
    #pconfig['ymax'] = float(T['y'].max() + 0.1*T['y'].max())
    #pconfig['ymin'] = float(T['y'].min() + 0.1*T['y'].min())

    data_labels = [{'name': 'UMAP1 vs UMAP2', 'xlab': 'UMAP1', 'ylab': 'UMAP2'}]
    if max_comp > 3:
        data_labels.append({'name': 'UMAP3 vs UMAP4', 'xlab': 'UMAP3', 'ylab': 'UMAP4'})
    pconfig['data_labels'] = data_labels

    data = []
    data1 = {}
    default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9',
                      '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']
    default_colors = sc.pl.palettes.default_20
    if color_by in df1.columns:
        groups = set(df1[color_by])
        if len(groups) > 10:
            default_colors = sc.pl.palettes.default_20
        elif len(groups) > 20:
            default_colors = sc.pl.palettes.default_28
        else:
            default_colors = default_102
        g_df = df1.groupby(color_by)
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
        if color_by in df2.columns:
            groups = set(df2[color_by])
            g_df = df2.groupby(color_by)
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
    parser = argparse.ArgumentParser(description="UMAP figure for QC report")
    parser.add_argument("exprs", help="scanpy h5ad file")
    parser.add_argument("--sample-info",
                        help="Optional sample sheet. Will subset expr table if needed",
                        dest="samples")
    parser.add_argument("--feature-info",
                        help="Optional feature info. Will subset expr table if needed",
                        dest="features")
    parser.add_argument("-o ",
                        "--output",
                        default="pca_mqc.png",
                        help="Output filename. Will default to pca_mqc.png, Optional [*_mqc.yaml]")
    parser.add_argument("--recipe", default="recipe_seurat",help="preprocessing recipe name as defined in Scanpy")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argparser()

    adata = sc.read(args.exprs)

    # standard preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_genes(adata, min_counts=20)
    sc.pp.normalize_total(adata, key_added='n_counts_all')
    f = sc.pp.filter_genes_dispersion(adata.X, flavor='cell_ranger', n_top_genes=1000, log=False)
    adata._inplace_subset_var(f.gene_subset)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    # neighbours
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata)   

    # cluster
    sc.tl.louvain(adata, resolution=0.8)

    # umap
    sc.tl.umap(adata)

    
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not E.var.columns.isin(S.index).all():
            raise ValueError("missing samples in sample info!")
        S = S.loc[E.var.columns, :]
    if args.features is not None:
        F = pd.read_csv(args.features, sep="\t", index_col=0)
        if not E.obs.index.isin(F.index).all():
            warnings.warn("missing annotations in feature info!")
        F = F.loc[E.obs.index, :]

    if args.output.endswith('_mqc.png'):
        pca_color = ['louvain', 'library_id']
        if 'Sample_Group' in adata.obs.columns:
            pca_color.append('Sample_Group')
        fig = sc.pl.umap(adata, return_fig=True, color=pca_color, ncols=1)
        fig.savefig(args.output, dpi=300, bbox_inches='tight')
    elif args.output.endswith('_mqc.yaml'):
        T = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names)
        multiqc_yaml(T, adata.obs, color_by='louvain')
    else:
        raise ValueError('Output file need to end with _mqc.png or _mqc.yaml')
