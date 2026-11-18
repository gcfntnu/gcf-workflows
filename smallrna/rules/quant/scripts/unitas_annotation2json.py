"""
Convert the annotation summary to json format (unitas.annotation_summary.txt)

    unitas_annotation2json.py SAMPLE/unitas.annotation_summary.txt --hits-per-target SAMPLE/unitas.hits_per_target.txt

misc_RNA is split by gene family (#379, unitas_misc_rna.py): {"misc_RNA": {"Y_RNA": ..., "7SK": ..., "7SL": ..., "other": ...}},
like the other classes with subclasses.
"""

import argparse
import sys
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

from unitas_misc_rna import check_total, misc_rna_families

def file2dict(fn):
    tree = {}
    with open(fn) as fh:
        txt = fh.read().splitlines()
    for line in txt:                                                      
        if not line.startswith('   '):
            p1, v = line.split('\t')
            tree[p1] = {}                    
        elif not line.startswith('      ') and line.startswith('   '):
            p2, v = line.split('\t')                            
            p2 = p2.strip()    
            tree[p1][p2] = {}                                            
        else:                      
            p3, v = line.split('\t')                            
            p3 = p3.strip()        
            tree[p1][p2][p3] = float(v)
    for line in txt:                                                      
        name, count = line.split('\t')
        name = name.strip()
        if name in tree and tree[name] == {}:
            tree[name] = float(count)
        else:
            for d, dd in tree.items():
                if isinstance(dd, dict):
                    if name in dd and dd[name] == {}:                         
                        tree[d][name] = float(count)
    return tree

def add_misc_rna_families(tree, hits):
    misc_rna = tree.get('misc_RNA', 0.0)
    if isinstance(misc_rna, dict):
        raise ValueError('{}: the annotation summary already has misc_RNA subclasses'.format(hits))
    counts = misc_rna_families(hits)
    check_total(counts, misc_rna, hits)
    if 'misc_RNA' in tree:
        tree['misc_RNA'] = counts
    return tree

if __name__ == '__main__':
    import json

    parser = argparse.ArgumentParser(description='Convert a UNITAS annotation summary to JSON')
    parser.add_argument('summary')
    parser.add_argument('--hits-per-target', required=True, dest='hits', help="the sample's unitas.hits_per_target.txt")
    args = parser.parse_args()
    try:
        tree = add_misc_rna_families(file2dict(args.summary), args.hits)
    except (OSError, ValueError) as error:
        sys.exit('unitas_annotation2json.py: {}'.format(error))
    json.dump(tree, sys.stdout)
