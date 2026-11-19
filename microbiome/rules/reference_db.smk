# -*- mode:snakemake -*-

ref = config['db']['reference_db']
DB_CONF = config['db'][ref]
if 'REF_DIR' not in locals() or not REF_DIR:
    REF_DIR = config.get('ref_dir')
if not REF_DIR:
    if ref == 'silva':
        include: join(GCFDB_DIR, 'silva.db')
    if ref == 'greegenes' or ref.lower() == 'gg':
        include: join(GCFDB_DIR, 'greengenes.db')
    if ref == 'unite':
        include: join(GCFDB_DIR, 'unite.db')

    REF_DIR = join(EXT_DIR, ref, DB_CONF['release'])

include:
    join(GCFDB_DIR, 'qiime2_classifiers.smk')
