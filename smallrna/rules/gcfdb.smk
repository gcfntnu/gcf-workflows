#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'reference_db.db')
include:
    join(GCFDB_DIR, 'spikein.db')
include:
    join(GCFDB_DIR, 'smallrna_calibrators.db')
include:
    join(GCFDB_DIR, 'contaminants.db')
include:
    join(GCFDB_DIR, 'rrna.db')
# unitas.db checks the organism; keep it before the other organism-specific references.
include:
    join(GCFDB_DIR, 'unitas.db')
include:
    join(GCFDB_DIR, 'mirbase.db')
