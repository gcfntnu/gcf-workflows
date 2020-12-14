#!/usr/bin/env python

import sys
import warnings
import glob
import os

import pandas as pd

if __name__ == '__main__':
    input_html = sys.argv[1]
    outdir = sys.argv[2]
    warnings.warn("INPUT {}".format(str(input)))
    filenames = glob.glob(os.path.join(os.path.split(input_html)[0], 'tables', '*.txt')) 
    for fn in filenames:
        warnings.warn("INPUT file: {}".format(fn))
        if fn.endswith('.html'):
            continue
        df = pd.read_csv(fn, sep='\t', index_col=0)
        pth, bn = os.path.split(fn)
        excel_bn = os.path.splitext(bn)[0] + '.xlsx'
        if not os.path.exists(os.path.join(pth, 'excel_format')):
            os.makedirs(os.path.join(pth, 'excel_format'))
        out = os.path.join(pth, 'excel_format', excel_bn)
        df.to_excel(out, index_label='ID', freeze_panes=[1,4])
        print('wrote file: {}'.format(out))
