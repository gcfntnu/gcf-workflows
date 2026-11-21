#-*- mode:snakemake -*-
"""mirflowz 0.10.0 as an optional quant_all route (quant.method: mirflowz, alone or combined, e.g. unitas,mirflowz).

mirflowz maps the shared filtered reads to the genome and transcriptome (segemehl, oligomap) and counts miRNAs,
isomiRs and pri-miRs against miRBase coordinates. It runs unmodified in gcfntnu/mirflowz (no inner containers), on the
reference prepared once by gcfdb/mirflowz.db. Delivery (bfq_all, multiqc_report) does not use it.

mirflowz always quality-filters and adapter-trims its input. On the already filtered reads that is made a no-op where
mirflowz allows it: reads passed as FASTA, so no quality filter, the adapter only trimmed where it occurs in full (overlap = adapter length, so
3' ends that match the start of the adapter are kept), the kit's own minimum length, and no N limit. mirflowz still
trims N bases at read ends (cutadapt --trim-n; it has no option to turn that off).
"""

MIRFLOWZ_SELECTED = 'mirflowz' in [m.strip() for m in str(config['quant']['method']).split(',')]

MIRFLOWZ_MEM_GB = int(config['quant']['mirflowz']['mem_gb'])
MIRFLOWZ_THREADS = int(config['quant']['mirflowz']['threads'])
if MIRFLOWZ_THREADS < 1:
    raise ValueError('quant.mirflowz.threads must be at least 1, got {}'.format(MIRFLOWZ_THREADS))

include:
    join(GCFDB_DIR, 'mirflowz.db')

MIRFLOWZ_MIN_LENGTH = filtered_min_length() or 1
MIRFLOWZ_OUT = join(QUANT_INTERIM, 'mirflowz')
# Keys of the run's mirflowz_config.yaml (written by run_mirflowz.py) that decide counting, recorded as they ran (#417)
MIRFLOWZ_RECORDED = ('minimum_length', 'overlap', 'max_n', 'max_length_reads', 'nh', 'extension', 'error_rate', 'q_value', 'p_value')


def mirflowz_settings():
    """The settings that decide mirflowz's counting besides MIRFLOWZ_RECORDED (quant_settings, #417)."""
    return [('adapter', config['adapter'], 'kit'),
            ('multimappers', 'divided between their loci by mirflowz (not configurable)', 'fixed')]

rule mirflowz_quant:
    input:
        fastq = expand(get_filtered_fastq(None)['R1'], sample=SAMPLES),
        prepared = rules.mirflowz_prepare.output.prepared,
        reference = rules.mirflowz_prepare.output.record,
        annotation = rules.mirflowz_annotation.output.annotation,
        annotation_record = rules.mirflowz_annotation.output.record
    output:
        mirna = join(MIRFLOWZ_OUT, 'all_mirna_counts.tab'),
        primir = join(MIRFLOWZ_OUT, 'all_pri-mir_counts.tab'),
        alignments = directory(join(MIRFLOWZ_OUT, 'alignments')),
        config = join(MIRFLOWZ_OUT, 'mirflowz_config.yaml')
    params:
        script = src_gcf('scripts/run_mirflowz.py'),
        samples = lambda wildcards, input: ' '.join('{}={}'.format(s, f) for s, f in zip(SAMPLES, input.fastq)),
        adapter = config['adapter'],
        min_length = MIRFLOWZ_MIN_LENGTH,
        out = MIRFLOWZ_OUT
    threads:
        MIRFLOWZ_THREADS
    resources:
        mem_mb = MIRFLOWZ_MEM_GB * 1024
    benchmark:
        join(QUANT_BENCHMARKS, 'mirflowz_quant.tsv')
    container:
        'docker://' + config['docker']['mirflowz']
    shell:
        'python {params.script} --prepared {input.prepared} --annotation {input.annotation} --adapter {params.adapter} '
        '--minimum-length {params.min_length} --threads {threads} --mem-mb {resources.mem_mb} --output {params.out} {params.samples}'

rule mirflowz_all:
    input:
        rules.mirflowz_quant.output.mirna,
        rules.mirflowz_quant.output.primir
