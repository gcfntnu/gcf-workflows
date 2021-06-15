import sys

import scanpy as sc
import scvelo as scv

A = sc.read(sys.argv[1])
B = sc.read(sys.argv[2])

merged = scv.utils.merge(A, B, copy=True)

merged.write(sys.argv[3])
