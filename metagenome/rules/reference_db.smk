# -*- mode:snakemake -*-

if 'REF_DIR' in locals() and REF_DIR is not None:
    pass
elif environ.get('GCF_REFDIR') is not None:
    REF_DIR = environ.get('GCF_REFDIR')
elif 'ref_dir' in config:
    REF_DIR = config['ref_dir']
elif 'reference_db' in config['db']:
    ref = config['db']['reference_db']
    if ref == 'langmead':
        include: join(GCFDB_DIR, 'langmead.db')
    elif ref == 'ncbi_16s':
        include: join(GCFDB_DIR, 'ncbi_16s.db')
    DB_CONF = config['db'][ref]
    REF_DIR = join(EXT_DIR, ref, DB_CONF['release'])
else:
    logger.error('Reference dir missing. Check config or other ways for ref dir setup')
    raise ValueError

if not 'ref' in locals():
    ref = os.path.dirname(REF_DIR)
    config['db']['reference_db'] = ref
    config['db'][ref] = {}

