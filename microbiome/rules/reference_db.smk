# -*- mode:snakemake -*-

if 'REF_DIR' in locals() and REF_DIR is not None:
    pass
elif environ.get('GCF_REFDIR') is not None:
    REF_DIR = environ.get('GCF_REFDIR')
elif 'ref_dir' in config:
    REF_DIR = config['ref_dir']
elif 'reference_db' in config['db']:
    ref = config['db']['reference_db']
    if ref == 'silva':
        include: join(GCFDB_DIR, 'silva.db')
    if ref == 'greegenes' or ref.lower() == 'gg':
        include: join(GCFDB_DIR, 'greengenes.db')
    if ref == 'unite':
        include: join(GCFDB_DIR, 'unite.db')

    DB_CONF = config['db'][ref]
    REF_DIR = join(EXT_DIR, ref, DB_CONF['release'])
else:
    logger.error('Reference dir missing. Check config or other ways for ref dir setup')
    raise ValueError

if not 'ref' in locals():
    ref = os.path.dirname(REF_DIR)
    config['db']['reference_db'] = ref
    config['db'][ref] = {}

#include:
#    join(GCFDB_DIR, 'silva.db')
#include:
#   join(GCFDB_DIR,  'greengenes.db')
include:
    join(GCFDB_DIR,  'qiime2_classifiers.smk')
