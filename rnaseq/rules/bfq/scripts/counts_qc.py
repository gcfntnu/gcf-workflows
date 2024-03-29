#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import pandas as pd
import numpy as np


def high_index(E, n_top=10):
    order = E.values.sum(1).argsort()[::-1]
    keep = order[:n_top]
    out = order[n_top:]
    genes = E.index[keep]
    other = E.index[out]
    return genes, other

def high_values(X, index, copy=True):
    x = X.loc[genes,:]
    if copy:
        x = X.copy()
    return X


def fig_expression(abundance, counts, feature_info, norm=True, msg=''):
    #abundance = 100* (abundance / abundance.sum(0))
    #counts = 100* (counts / counts.sum(0))
    a_df = pd.concat([abundance, feature_info], axis=1)
    c_df = pd.concat([counts, feature_info], axis=1)
    keys = {}
    default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c',
                      '#8085e9','#f15c80', '#e4d354', '#2b908f',
                      '#f45b5b', '#91e8e1', '#bdbdbd', '#434348', '#90ed7d', '#f7a35c',
                      '#8085e9','#f15c80', '#e4d354', '#2b908f',
                      '#f45b5b', '#91e8e1', '#bdbdbd']
    for i, k in enumerate(a_df.index):
        gene_name = a_df.loc[k, 'gene_name']
        if pd.isnull(gene_name):
            gene_name = k
        gene_biotype =  str(a_df.loc[k, 'gene_biotype'])
        if gene_biotype == 'other':
            name = 'Other'
            color = '#bdbdbd'
        else:
            gene_biotype = gene_biotype.replace('_', ' ').title()
            name = str(gene_name) + ' (' + gene_biotype + ')' 
            color = default_colors[i]
        name = name.replace('_', ' ').title()
        keys[k] = {'color': color, 'name': name}
    
    
    # Config for the plot
    pconfig = {
        'id': 'gene_high',
        'title': 'Top 10 most abundant genes',
        'ylab': '# Reads',
        'tt_decimals': 1,
        'tt_percentages': False,
        'cpswitch': False,
        'data_labels' : [{'name': 'counts', 'ylab': '# Reads'}, {'name': 'TPM', 'ylab': 'TPM fraction'}]
    }
    
    section = {}
    section['id'] = 'gene_high'
    section['section_name'] = 'Abundant Genes'
    section['description'] = '. Sorted by average TPM expression across all samples.' + msg
    section['plot_type'] = 'bargraph'
    
    section['pconfig'] = pconfig
    section['categories'] = [keys, keys]

    data = []
    a_dict = {}
    for k, v in a_df.to_dict().items():
        if not k in ['gene_name', 'gene_biotype']:
            a_dict[k] = v
    c_dict = {}
    for k, v in c_df.to_dict().items():
        if not k in ['gene_name', 'gene_biotype']:
            c_dict[k] = v
    
    section['data'] = [c_dict, a_dict]

    return section


def fig_biotype(E, C, F, rel=0.01, msg=''):
    def setup_df(X, F, rel=0.01, frac=0.5):
        if frac > 1 :
            n_samples_lim = int(frac)
        else:
            n_samples_lim = np.floor(E.shape[1] * float(frac))
        df = pd.concat([X, F["gene_biotype"]], axis=1)
        tab = df.groupby("gene_biotype").sum()
        tab_rel = tab / tab.sum(0)
        keep = (tab_rel >= rel).sum(1) > n_samples_lim
        if keep.sum() == 0:
            #fixme
            raise ValueError
        other = tab.loc[~keep,:].sum(0)
        other = pd.DataFrame(other).T
        other.index = ['other']
        tab = tab.loc[keep, :]
        tab = tab.append(other)
        order = tab.sum(1).argsort()[::-1]
        tab = tab.iloc[order, :]
        return tab
    
    a_tab = setup_df(E, F, rel)
    c_tab = setup_df(C, F, rel)
    cats = []
    keys = {}
    default_colors = ["#7cb5ec", "#434348", "#90ed7d", "#f7a35c", "#8085e9", "#f15c80",
                      "#e4d354", "#2b908f", "#f45b5b", "#91e8e1"]
    for i, k in enumerate(a_tab.index):
        if k == k.upper():
            name = k
        else:
            name = k.replace("_", " ").title()
        if k  == 'other':
            color = "#bdbdbd"
        else:
            color = default_colors[i]
        keys[k] = {"color": color, "name": name}
    # Config for the plot
    pconfig = {
        "id": "gene_biotypes",
        "title": "Gene Biotype Counts",
        'ylab': '# Reads',
        'tt_decimals': 1,
        'tt_percentages': False,
        'cpswitch': False,
        'data_labels' : [{'name': 'counts', 'ylab': '# Reads'}, {'name': 'TPM', 'ylab': 'TPM'}]
    }

    section = {}
    section["id"] = "gene_biotypes"
    section["section_name"] = "Gene Biotypes Count"
    section["description"] = " .Summary of gene biotype annotation. Low abundant biotypes (<1%) is put into `other` category. {}".format(msg)
    section["plot_type"] = "bargraph"

    section["pconfig"] = pconfig
    section["categories"] = [keys, keys]
    section["data"] = [c_tab.to_dict(), a_tab.to_dict()]

    return section


def argparser():
    parser = argparse.ArgumentParser(description="Higly expressed genes figure for QC report")
    parser.add_argument("--figure", '-f', help="Figure type ", dest='fig', choices=['top_genes', 'top_biotypes'])
    parser.add_argument("--abundance", help="TPM values", dest='exprs')
    parser.add_argument("--counts", help="Count values", dest='counts')
    parser.add_argument("--sample-info", help="Optional sample sheet. Will subset expr table if needed",
                        dest="samples")
    parser.add_argument("--feature-info", help="Required feature info. Will subset expr table if needed",
                        dest="features", required=True)
    parser.add_argument("--gene-lengths", help="Optional gene length table. See `--gene-length-cutoff`",
                        dest="lengths", required=False)
    parser.add_argument("--gene-length-cutoff", type=int, default=100, dest='len_cutoff', 
                        help="Discard genes with effective gene length below this value in > 50% of samples")
    parser.add_argument("-o ", "--output", default="biotypes_mqc.yaml",
                        help="Output filename. Will default to biotypes_mqc.yaml")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    msg = ''
    args = argparser()
    E = pd.read_csv(args.exprs, sep="\t", index_col=0)
    E.columns = keep_samples = E.columns.astype(str)
    #E = E / 1E4 # normalize TPM to 0-100
    keep_features = E.index

    C = pd.read_csv(args.counts, sep="\t", index_col=0)
    C.columns = C.columns.astype(str)

    #C = 100*(C / C.sum(0))
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not E.columns.isin(S.index).all():
            raise ValueError("missing samples in sample info!")
        keep_samples = S.index

    E = E.loc[:,keep_samples]
    C = C.loc[:,keep_samples]
    
    if args.lengths is not None:
        L = pd.read_csv(args.lengths, sep="\t", index_col=0)
        L = L.loc[:,keep_samples]
        if not keep_features.isin(L.index).all():
            raise ValueError("missing gene length for some genes!")
        L = L.loc[keep_features,:]
        if args.len_cutoff is not None and args.len_cutoff > 0:
            n_samples = E.shape[1]
            out = (L <= args.len_cutoff).sum(1) > (n_samples / 2.)
            if sum(out) > 0:
                keep_features = keep_features[out == False]
                msg += '. {} features removed from calculation due to low effective gene length (<{})'.format(sum(out), args.len_cutoff)
                print(msg)
    
    E = E.loc[keep_features,:]
    C = C.loc[keep_features,:]
    L = L.loc[keep_features,:]
    
    F = pd.read_csv(args.features, sep="\t", index_col='gene_id')
    if not keep_features.isin(F.index).all():
        warnings.warn("missing annotations in feature info!")
        keep_features = list(set(keep_features).intersection(F.index))
        E = E.loc[keep_features,:]
        C = C.loc[keep_features,:]
    F = F.loc[keep_features, :]
    if not 'gene_biotype' in F.columns:
        if 'gene_type' in F.columns: #gencode
            F = F.rename(columns={'gene_type': 'gene_biotype'})
        else:
            print(F.head(n=2))
            raise ValueError('Feature info needs column `gene_biotype`')   
    if 'gene_name' not in F.columns:
        if 'gene_id' in F.columns:
            F['gene_name'] = F['gene_id'].copy()
        else:
            F['gene_name'] = F.index.values.copy()
    F = F[['gene_name', 'gene_biotype']]

    # renormalize tpm abundance
    A = C / L
    E = (A/A.sum(0)) 
    keep, inv_keep = high_index(E, 10)

    abundance = E.loc[keep,:]
    other = E.loc[inv_keep,:].sum(0)
    other.name = "other"
    abundance = abundance.append(other)
    counts = C.loc[keep,:]
    other = C.loc[inv_keep,:].sum(0)
    other.name = "other"
    counts = counts.append(other)
    feature_info = F.loc[keep,:]
    other = pd.Series(['other','other'], name="other", index=feature_info.columns)
    feature_info = feature_info.append(other)

    
    if args.fig == 'top_genes':
        section =  fig_expression(abundance, counts, feature_info, msg=msg)
    elif args.fig == 'top_biotypes':
        section =  fig_biotype(E, C, F, msg=msg)
    elif args.fig == 'top_influencers':
        pass
    elif args.fig == 'top_misfits':
        pass
    else:
        print(args)
        raise ValueError
    
    
    with open(args.output, 'w') as fh:
        yaml.dump(section, fh, default_flow_style=False, sort_keys=False)
