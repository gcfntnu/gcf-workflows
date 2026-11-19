# -*- mode:snakemake -*-

include: join(GCFDB_DIR, 'kraken_settings.smk')

provider = config['db']['reference_db']
if provider not in ('langmead', 'ncbi_16s'):
    raise ValueError(f"Unsupported metagenome db.reference_db: {provider!r}")
DB_CONF = config['db'][provider]
BRACKEN_N_MERS = bracken_read_lengths(provider, DB_CONF['assembly'])

if PROJECT_KRAKEN_DB:
    REF_DIR = PROJECT_KRAKEN_DB
    KRAKEN_INDEX = kraken_index_files(REF_DIR)
elif provider == 'langmead':
    include: join(GCFDB_DIR, 'langmead.db')
    REF_DIR = dirname(rules.langmead_kraken_prebuild.output.hash)
    KRAKEN_INDEX = rules.langmead_shmem.output
else:
    include: join(GCFDB_DIR, 'ncbi_16s.db')
    REF_DIR = dirname(rules.ncbi_16s_kraken_build.output.hash)
    KRAKEN_INDEX = rules.ncbi_16s_shmem.output

def metagenome_database_input(path):
    # An explicit directory is operator-supplied, even if it happens to match a getter's output path.
    if PROJECT_KRAKEN_DB and not os.path.isfile(path):
        raise ValueError(f"Missing required database file in project ref_dir: {path}")
    return path
