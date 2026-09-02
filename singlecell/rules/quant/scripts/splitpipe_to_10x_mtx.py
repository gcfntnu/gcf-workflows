#!/usr/bin/env python

import argparse
import os
import shutil
import gzip

import pandas as pd
import scipy
import scipy.io

def argparser():
    parser = argparse.ArgumentParser(description='Parse mtx conversion to 10xgenomics type (v2) mtx ')
    parser.add_argument('-i', '--input', required=True, help="parse mtx file")
    parser.add_argument('--barcode-info', help='Optional barcode info.', dest='barcodes', default=None)
    parser.add_argument('--feature-info', help='Optional feature info.', dest='features', default=None)
    parser.add_argument('--compress', action='store_true', help='gzip compress output files (.gz)')
    parser.add_argument('--enforce-float', action='store_true', help='store integer counts as floats (f32)')
    parser.add_argument('-o ', '--output', required=True, help='Output filename (matrix.mtx.gz)')
    args = parser.parse_args()
    return args


def identify_feature_info(pth):
    accepted = set("all_genes.csv target_genes.csv all_guides.csv".split())
    for fn in os.listdir(pth):
        if os.path.basename(fn) in accepted:
            return os.path.basename(fn)
    raise FileNotFoundError
        
if __name__ == '__main__':
    args = argparser()
    compress = os.path.splitext(args.output)[-1] == '.gz'
    args.compress = (args.compress) or compress
    mtx_dirname = os.path.dirname(args.input)
    mtx = None
    mtx = scipy.io.mmread(args.input).T
    if args.enforce_float:
        mtx = mtx.asfptype()

    if args.barcodes is None:
        args.barcodes = os.path.join(mtx_dirname, "cell_metadata.csv")
    if args.features is None:
        fn = identify_feature_info(mtx_dirname)
        args.features = os.path.join(mtx_dirname, fn) 
        
    barcodes = pd.read_csv(args.barcodes)
    barcodes = barcodes.loc[:,"bc_wells"]
    features = pd.read_csv(args.features)
    features = features[['gene_id', 'gene_name']]
    if mtx is None:
        # copy mtx file
        assert(args.input != args.output)
        shutil.copyfile(args.input, args.output)
    else:
        n_rows, n_cols = mtx.shape
        parse_comment = f"Rows=genes ({n_rows}), Cols=cells ({n_cols}) "
        if args.compress: 
            mtx_fh = gzip.open(args.output, 'wb')
            scipy.io.mmwrite(mtx_fh, mtx, comment=parse_comment, precision=2)
        else:
            scipy.io.mmwrite(args.output, mtx, comment=parse_comment, precision=2)

    output_dirname = os.path.dirname(args.output)
    if args.compress:
        barcodes_output_fn = os.path.join(output_dirname, "barcodes.tsv.gz")
        features_output_fn = os.path.join(output_dirname, "genes.tsv.gz")
        barcodes.to_csv(barcodes_output_fn, compression='gzip', header=False, index=False)
        features.to_csv(features_output_fn, compression='gzip', header=False, index=False, sep="\t")
    else:
        barcodes_output_fn = os.path.join(output_dirname, "barcodes.tsv")
        features_output_fn = os.path.join(output_dirname, "genes.tsv")        
        barcodes.to_csv(barcodes_output_fn, header=False, index=False)
        features.to_csv(features_output_fn, header=False, index=False, sep="\t")
