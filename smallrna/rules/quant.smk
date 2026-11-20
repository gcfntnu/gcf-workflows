#-*- mode:snakemake -*-

QUANT_METHODS = ('unitas', 'mirge3', 'mirflowz', 'tesmall', 'tinyrna', 'smrnaseq', 'miraligner')
QUANT_SELECTED = [m.strip() for m in str(config['quant']['method']).split(',')]
# Set before the routes are included: the TEsmall database and the tiny-count reference are built on the mirflowz
# reference (gcfdb/mirflowz.db), and tiny-count uses TEsmall's genome index (gcfdb/tesmall.db).
TESMALL_SELECTED = 'tesmall' in QUANT_SELECTED
TINYRNA_SELECTED = 'tinyrna' in QUANT_SELECTED
# Snakemake benchmark files of the routes' rules (runtime and memory, read by the benchmark's evaluate, #386):
# <rule>/<sample>.tsv per sample, <rule>.tsv per project. The shared reference builds of gcfdb/mirflowz.db and
# gcfdb/tesmall.db write theirs in the reference they built (<reference>/benchmarks/<rule>.tsv), not here.
QUANT_BENCHMARKS = join(QUANT_INTERIM, 'benchmarks')

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
include:
    'quant/miraligner.smk'

RECORDED_ROUTES = ('mirflowz', 'tesmall', 'tinyrna', 'smrnaseq', 'miraligner')


def quant_reference_record(route):
    """(miRNA database, release, source, reference) a route annotated with: the mirflowz reference's for mirflowz, TEsmall
    and tiny-count (db.mirbase.release, or NA for a db.mirflowz.custom reference, whose annotation is the project's own),
    gcfdb's miRBase files for nf-core/smrnaseq, and the database written from them for miraligner."""
    if route == 'smrnaseq':
        return 'mirbase', MIRBASE_RELEASE, 'db.mirbase.release', MIRBASE_EXT
    if route == 'miraligner':
        return 'mirbase', MIRBASE_RELEASE, 'db.mirbase.release', MIRALIGNER_DB
    if MIRFLOWZ_CUSTOM:
        return 'NA', 'NA', 'db.mirflowz.custom', MIRFLOWZ_REF
    return 'mirbase', MIRFLOWZ_MIRBASE, 'db.mirbase.release', MIRFLOWZ_REF


rule quant_reference:
    """Which miRBase release a route annotated with (quant/<route>/reference.tsv, the layout of UNITAS's and miRge3's
    records), read by the benchmark before scoring (#407). The route's count table is its input, so the record is
    written only after that table was (re)built for the current configuration, and it carries the table's SHA-256
    (table_sha256), so the benchmark accepts it only for that very table."""
    input:
        lambda wildcards: {'mirflowz': rules.mirflowz_quant.output.mirna, 'tesmall': rules.tesmall_counts.output[0],
                           'tinyrna': rules.tinyrna_counts.output[0], 'smrnaseq': rules.smrnaseq_counts.output.mirna,
                           'miraligner': rules.miraligner_counts.output[0]}[wildcards.route]
    output:
        join(QUANT_INTERIM, '{route}', 'reference.tsv')
    params:
        database = lambda wildcards: quant_reference_record(wildcards.route)[0],
        release = lambda wildcards: quant_reference_record(wildcards.route)[1],
        source = lambda wildcards: quant_reference_record(wildcards.route)[2],
        reference = lambda wildcards: quant_reference_record(wildcards.route)[3]
    wildcard_constraints:
        route = '|'.join(RECORDED_ROUTES)
    container:
        'docker://' + config['docker']['default']
    shell:
        "sha=$(sha256sum {input:q} | cut -d ' ' -f 1) && "
        "printf 'field\\tvalue\\tsource\\nmirna_database\\t%s\\t%s\\nmirna_release\\t%s\\t%s\\nreference\\t%s\\tconfiguration\\ntable_sha256\\t%s\\t%s\\n' "
        "{params.database:q} {params.source:q} {params.release:q} {params.source:q} {params.reference:q} \"$sha\" {input:q} > {output:q}"


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
        elif m == 'miraligner':
            out.extend(rules.miraligner_all.input)
    out.extend(join(QUANT_INTERIM, m, 'reference.tsv') for m in methods if m in RECORDED_ROUTES)
    return out

rule quant_all:
    input:
        quant_all
