"""merge scanpy objects


"""
import sys

import numpy as np
import scanpy as sc
import scvelo as scv

A = sc.read(sys.argv[1])
B = sc.read(sys.argv[2])

merged = scv.utils.merge(A, B, copy=True)

# add nuclear fraction to obs data
exon_sum = merged.layers['spliced'].sum(axis=1)
intron_sum = merged.layers['unspliced'].sum(axis=1)
nuclear_fraction = intron_sum/(exon_sum + intron_sum)
if hasattr(nuclear_fraction, "A1"):
    nuclear_fraction = nuclear_fraction.A1
merged.obs['nuclear_fraction'] = nuclear_fraction


merged.write(sys.argv[3])
