#-*- mode:snakemake -*-
import json

QUANT_METHODS = ('unitas', 'mirge3', 'mirflowz', 'tesmall', 'tinyrna', 'smrnaseq', 'biocore_htseq', 'biocore_miraligner')
QUANT_SELECTED = [m.strip() for m in str(config['quant']['method']).split(',')]
# Set before the routes are included: the TEsmall database and the tiny-count and biocore_htseq references are built on the
# mirflowz reference (gcfdb/mirflowz.db), and tiny-count and biocore_htseq use TEsmall's genome (gcfdb/tesmall.db).
TESMALL_SELECTED = 'tesmall' in QUANT_SELECTED
TINYRNA_SELECTED = 'tinyrna' in QUANT_SELECTED
BIOCORE_HTSEQ_SELECTED = 'biocore_htseq' in QUANT_SELECTED
# Snakemake benchmark files of the routes' rules (runtime and memory, read by the benchmark's evaluate, #386):
# <rule>/<sample>.tsv per sample, <rule>.tsv per project. The shared reference builds of gcfdb/mirflowz.db,
# gcfdb/tesmall.db, gcfdb/biocore_htseq.db and gcfdb/biocore_miraligner.db write theirs in the reference they built
# (<reference>/benchmarks/<rule>.tsv), not here.
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
    'quant/biocore_htseq.smk'
include:
    'quant/biocore_miraligner.smk'

RECORDED_ROUTES = ('mirflowz', 'tesmall', 'tinyrna', 'smrnaseq', 'biocore_htseq', 'biocore_miraligner')


def quant_reference_record(route):
    """(miRNA database, release, source, reference, GFF identity) a route annotated with: the mirflowz reference's for mirflowz,
    TEsmall, tiny-count and biocore_htseq (db.mirbase.release, or NA for a db.mirflowz.custom reference, whose annotation is the
    project's own), gcfdb's miRBase files for nf-core/smrnaseq, and the database written from them for biocore_miraligner. The
    GFF identity names the files built from the miRNA GFF (#440): the corrected miRBase GFF3's, or `custom`."""
    if route == 'smrnaseq':
        return 'mirbase', MIRBASE_RELEASE, 'db.mirbase.release', MIRBASE_EXT, MIRBASE_GFF_ID
    if route == 'biocore_miraligner':
        return 'mirbase', MIRBASE_RELEASE, 'db.mirbase.release', BIOCORE_MIRALIGNER_DB, MIRBASE_GFF_ID
    if MIRFLOWZ_CUSTOM:
        return 'NA', 'NA', 'db.mirflowz.custom', MIRFLOWZ_REF, MIRFLOWZ_GFF_ID
    return 'mirbase', MIRFLOWZ_MIRBASE, 'db.mirbase.release', MIRFLOWZ_REF, MIRFLOWZ_GFF_ID


rule quant_reference:
    """Which miRBase release a route annotated with (quant/<route>/reference.tsv, the layout of UNITAS's and miRge3's
    records), read by the benchmark before scoring (#407). The route's count table is its input, so the record is
    written only after that table was (re)built for the current configuration, and it carries the table's SHA-256
    (table_sha256), so the benchmark accepts it only for that very table."""
    input:
        lambda wildcards: {'mirflowz': rules.mirflowz_quant.output.mirna, 'tesmall': rules.tesmall_counts.output[0],
                           'tinyrna': rules.tinyrna_counts.output[0], 'smrnaseq': rules.smrnaseq_counts.output.mirna,
                           'biocore_htseq': rules.biocore_htseq_counts.output[0],
                           'biocore_miraligner': rules.biocore_miraligner_counts.output[0]}[wildcards.route]
    output:
        join(QUANT_INTERIM, '{route}', 'reference.tsv')
    params:
        database = lambda wildcards: quant_reference_record(wildcards.route)[0],
        release = lambda wildcards: quant_reference_record(wildcards.route)[1],
        source = lambda wildcards: quant_reference_record(wildcards.route)[2],
        reference = lambda wildcards: quant_reference_record(wildcards.route)[3],
        gff_id = lambda wildcards: quant_reference_record(wildcards.route)[4]
    wildcard_constraints:
        route = '|'.join(RECORDED_ROUTES)
    container:
        'docker://' + config['docker']['default']
    shell:
        "sha=$(sha256sum {input:q} | cut -d ' ' -f 1) && "
        "printf 'field\\tvalue\\tsource\\nmirna_database\\t%s\\t%s\\nmirna_release\\t%s\\t%s\\nreference\\t%s\\tconfiguration\\n"
        "mirna_gff_id\\t%s\\tconfiguration\\ntable_sha256\\t%s\\t%s\\n' "
        "{params.database:q} {params.source:q} {params.release:q} {params.source:q} {params.reference:q} {params.gff_id:q} \"$sha\" {input:q} "
        "> {output:q}"


def quant_table(route):
    """A route's count table, which the benchmark reads (adapters.TABLES)."""
    return {'unitas': rules.unitas_mirtable.output.counts, 'mirge3': rules.mirge3_counts.output[0], 'mirflowz': rules.mirflowz_quant.output.mirna,
            'tesmall': rules.tesmall_counts.output[0], 'tinyrna': rules.tinyrna_counts.output[0], 'smrnaseq': rules.smrnaseq_counts.output.mirna,
            'biocore_htseq': rules.biocore_htseq_counts.output[0], 'biocore_miraligner': rules.biocore_miraligner_counts.output[0]}[route]


ROUTE_SETTINGS = {'unitas': unitas_settings, 'mirge3': mirge3_settings, 'mirflowz': mirflowz_settings, 'tesmall': tesmall_settings,
                  'tinyrna': tinyrna_settings, 'smrnaseq': smrnaseq_settings, 'biocore_htseq': biocore_htseq_settings,
                  'biocore_miraligner': biocore_miraligner_settings}


rule quant_settings:
    """The settings that decide a route's counting (quant/<route>/settings.tsv, #417), read by the benchmark with its
    results: each route's *_settings() in rules/quant/<route>.smk, built from the variables its commands use, and for
    mirflowz the keys MIRFLOWZ_RECORDED of the configuration it ran with. Written after the count table, with its SHA-256."""
    input:
        table = lambda wildcards: quant_table(wildcards.route),
        config = lambda wildcards: rules.mirflowz_quant.output.config if wildcards.route == 'mirflowz' else []
    output:
        join(QUANT_INTERIM, '{route}', 'settings.tsv')
    params:
        script = src_gcf('quant/scripts/quant_settings.py'),
        settings = lambda wildcards: json.dumps([list(row) for row in ROUTE_SETTINGS[wildcards.route]()]),
        recorded = lambda wildcards, input: '--recorded {}:{}'.format(input.config, ','.join(MIRFLOWZ_RECORDED)) if input.config else ''
    wildcard_constraints:
        route = '|'.join(QUANT_METHODS)
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} --table {input.table:q} --settings {params.settings:q} {params.recorded} --output {output:q}'


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
        elif m == 'biocore_htseq':
            out.extend(rules.biocore_htseq_all.input)
        elif m == 'biocore_miraligner':
            out.extend(rules.biocore_miraligner_all.input)
    out.extend(join(QUANT_INTERIM, m, 'reference.tsv') for m in methods if m in RECORDED_ROUTES)
    out.extend(join(QUANT_INTERIM, m, 'settings.tsv') for m in methods)  # each route's counting settings (#417)
    return out

rule quant_all:
    input:
        quant_all
