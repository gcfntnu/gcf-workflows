import sys
import pyensembl

docker = sys.argv[1]

line = 'pyensembl,{},{}'.format(pyensembl.__version__, docker)
sys.stdout.write(line)
