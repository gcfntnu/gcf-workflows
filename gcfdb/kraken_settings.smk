# -*- mode:snakemake -*-
from os.path import realpath

# Bind project database aliases while planning. Kraken2 and Bracken can run at different times.
PROJECT_KRAKEN_DB = realpath(config['ref_dir']) if config.get('ref_dir') else None

def langmead_database_directory(root):
    settings = config['db']['langmead']
    return join(root, 'langmead', f"release-{settings['release']}", 'metagenome', settings['assembly'] or '')

def bracken_read_lengths(provider, assembly):
    if provider == 'langmead' and assembly.startswith('k2'):
        return [50, 75, 100, 150, 200, 250, 300]
    return [100, 150, 200, 250]

def kraken_index_files(directory):
    return [join(directory, name) for name in ('hash.k2d', 'opts.k2d', 'taxo.k2d')]
