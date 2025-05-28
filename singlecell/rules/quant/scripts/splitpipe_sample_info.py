import sys
import yaml
import pandas as pd


config_fn, cell_metadata_unfiltered_fn = sys.argv[1:]
with open("config.yaml") as fh:
    conf = yaml.safe_load(fh)
sample_info = pd.DataFrame.from_dict(conf['wells'], orient='index')
cell_metadata_unfiltered = pd.read_csv(cell_metadata_unfiltered_fn)[['sample', 'bc_wells']]
cell_metadata_unfiltered.columns = ['Sample_ID', 'barcode']
out = cell_metadata_unfiltered.merge(sample_info, on='Sample_ID', how='left').set_index('barcode')
out.to_csv(sys.stdout, sep="\t")
    
