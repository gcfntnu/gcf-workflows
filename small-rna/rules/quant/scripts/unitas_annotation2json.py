"""
Convert the annotation summary to json format (unitas.annotation_summary.txt)
"""

import sys
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

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

if __name__ == '__main__':
    import json

    tree = file2dict(sys.argv[1])
    json.dump(tree, sys.stdout)
