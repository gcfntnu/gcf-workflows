# -*- mode:snakemake -*-

provider = config['db']['reference_db']
if provider == 'langmead':
    include: join(GCFDB_DIR, 'langmead.db')
    REF_DIR = dirname(rules.langmead_kraken_prebuild.output.hash)
    KRAKEN_INDEX = rules.langmead_shmem.output
    BRACKEN_N_MERS = LM_BRACKEN_N_MERS
elif provider == 'ncbi_16s':
    include: join(GCFDB_DIR, 'ncbi_16s.db')
    REF_DIR = dirname(rules.ncbi_16s_kraken_build.output.hash)
    KRAKEN_INDEX = rules.ncbi_16s_shmem.output
    BRACKEN_N_MERS = NCBI_BRACKEN_N_MERS
else:
    raise ValueError(f"Unsupported metagenome db.reference_db: {provider!r}")

DB_CONF = config['db'][provider]

if config.get('ref_dir'):
    REF_DIR = config['ref_dir']
    KRAKEN_INDEX = [join(REF_DIR, name) for name in ('hash.k2d', 'opts.k2d', 'taxo.k2d')]

def metagenome_database_input(path):
    # An explicit directory is operator-supplied, even if it happens to match a getter's output path.
    if config.get('ref_dir') and not os.path.isfile(path):
        raise ValueError(f"Missing required database file in project ref_dir: {path}")
    return path
