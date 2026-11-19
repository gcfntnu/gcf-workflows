#-*- mode:snakemake -*-

QUANT_METHODS = ('unitas', 'mirge3', 'mirflowz', 'tesmall', 'tinyrna', 'smrnaseq')
QUANT_SELECTED = [m.strip() for m in str(config['quant']['method']).split(',')]
# Set before the routes are included: the TEsmall database and the tiny-count reference are built on the mirflowz
# reference (gcfdb/mirflowz.db), and tiny-count uses TEsmall's genome index (gcfdb/tesmall.db).
TESMALL_SELECTED = 'tesmall' in QUANT_SELECTED
TINYRNA_SELECTED = 'tinyrna' in QUANT_SELECTED

include:
    'quant/mirge3.smk'
include:
    'quant/unitas.smk'
include:
    'quant/mirflowz.smk'
include:
    'quant/tesmall.smk'
include:
    'quant/tinyrna.smk'
include:
    'quant/smrnaseq.smk'

GENOME_ROUTES = ('mirflowz', 'tesmall', 'tinyrna')


rule quant_reference:
    """Which miRBase release a route on the mirflowz reference annotated with (quant/<route>/reference.tsv, the layout of
    UNITAS's and miRge3's records), read by the benchmark before scoring (#407): db.mirbase.release, or NA for a
    db.mirflowz.custom reference, whose annotation is the project's own. The route's count table is its input, so the
    record is written only after that table was (re)built for the current configuration and never stamps an older one."""
    input:
        lambda wildcards: {'mirflowz': rules.mirflowz_quant.output.mirna, 'tesmall': rules.tesmall_counts.output[0],
                           'tinyrna': rules.tinyrna_counts.output[0]}[wildcards.route]
    output:
        join(QUANT_INTERIM, '{route}', 'reference.tsv')
    params:
        database = 'NA' if MIRFLOWZ_CUSTOM else 'mirbase',
        release = 'NA' if MIRFLOWZ_CUSTOM else MIRFLOWZ_MIRBASE,
        source = 'db.mirflowz.custom' if MIRFLOWZ_CUSTOM else 'db.mirbase.release',
        reference = MIRFLOWZ_REF
    wildcard_constraints:
        route = '|'.join(GENOME_ROUTES)
    container:
        'docker://' + config['docker']['default']
    shell:
        "printf 'field\\tvalue\\tsource\\nmirna_database\\t%s\\t%s\\nmirna_release\\t%s\\t%s\\nreference\\t%s\\tconfiguration\\n' "
        "{params.database:q} {params.source:q} {params.release:q} {params.source:q} {params.reference:q} > {output:q}"


def quant_all(wildcards):
    methods = [m.strip() for m in str(config['quant']['method']).split(',')]
    unknown = [m for m in methods if m not in QUANT_METHODS]
    if unknown:
        raise ValueError('quant.method: unsupported value(s) {}; choose from {}'.format(', '.join(unknown), ', '.join(QUANT_METHODS)))
    out = []
    for m in methods:
        if m == 'unitas':
            out.extend(rules.unitas_all.input)
        elif m == 'mirge3':
            out.extend(rules.mirge3_all.input)
        elif m == 'mirflowz':
            out.extend(rules.mirflowz_all.input)
        elif m == 'tesmall':
            out.extend(rules.tesmall_all.input)
        elif m == 'tinyrna':
            out.extend(rules.tinyrna_all.input)
        elif m == 'smrnaseq':
            out.extend(rules.smrnaseq_all.input)
    out.extend(join(QUANT_INTERIM, m, 'reference.tsv') for m in methods if m in GENOME_ROUTES)
    return out

rule quant_all:
    input:
        quant_all
