import sys

import csv
import os
import subprocess


tab = csv.reader(open(sys.argv[1]))
header = next(tab)
bam_files = []
for row in tab:
    name, mol_h5 = row
    pth = os.path.dirname(mol_h5)
    bam_fn = os.path.join(pth, 'possorted_genome_bam.bam')
    assert(os.path.exists(bam_fn))
    bam_files.append(bam_fn)
cmd = ['sambamba', 'merge', '-t' ,'8', sys.argv[2]]
cmd.extend(bam_files)
print(' '.join(cmd))
subprocess.run(cmd)
