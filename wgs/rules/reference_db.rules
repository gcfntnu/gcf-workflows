if 'REF_DIR' in locals() and REF_DIR is not None:
    pass
elif environ.get('GCF_REFDIR') is not None:
    REF_DIR = environ.get('GCF_REFDIR')
elif 'ref_dir' in config:
    REF_DIR = config['ref_dir']
elif 'reference_db' in config['db']:
    ref = config['db']['reference_db']
    DB_CONF = config['db'][ref]
    #REF_DIR = join(EXT_DIR, ref, DB_CONF.get('release', ''))
    REF_DIR = join(EXT_DIR, ref)
else:
    logger.error('Reference dir missing. Check config or other ways for ref dir setup')
    raise ValueError

include:
    join(GCFDB_DIR, 'indexes.rules')
include:
    join(GCFDB_DIR, 'ncbi.db')
include:
    join(GCFDB_DIR, 'opengene.db')

