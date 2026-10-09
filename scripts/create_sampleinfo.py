import sys

import pandas as pd
import peppy
import yaml


def col_is_empty(col):
    if all(col == ''):
        return True
    if col.isna().all():
        return True
    if col.isnull().all():
        return True
    return False


def sampleinfo_from_yaml(config):
    samples = config.get('samples')
    if not isinstance(samples, dict):
        raise ValueError("Native workflow config must contain a top-level 'samples' mapping")

    df = pd.DataFrame.from_dict(samples, orient='index')
    df.index.name = 'Sample_ID'
    if 'Sample_ID' not in df.columns:
        df = df.reset_index()

    cols = list(df.columns)
    cols.remove('Sample_ID')
    cols.insert(0, 'Sample_ID')
    df = df.loc[:, cols]

    empty = df.apply(col_is_empty, axis=0)
    return df.loc[:, ~empty]


def sampleinfo_from_peppy(fn):
    pep = peppy.Project(fn)
    df = pep.sample_table
    df = df.loc[:, df.convert_dtypes().dtypes != 'object']
    if 'sample_name' in df.columns and 'Sample_ID' not in df.columns:
        df = df.rename(columns={'sample_name': 'Sample_ID'})
    return df


def sampleinfo(fn):
    with open(fn) as fh:
        config = yaml.safe_load(fh)

    if isinstance(config, dict) and 'samples' in config:
        return sampleinfo_from_yaml(config)

    return sampleinfo_from_peppy(fn)


if __name__ == '__main__':
    fn = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else sys.stdout
    sampleinfo(fn).to_csv(out, sep='\t', index=False)
