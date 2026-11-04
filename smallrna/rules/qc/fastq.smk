#-*- mode: snakemake -*-
"""miRTrace quality control on the merged, untrimmed FASTQs.

miRTrace trims adapters itself, using the configured adapter and the kit's protocol.
"""

MIRTRACE_DIR = join(QC_INTERIM, 'mirtrace')
MIRTRACE_SPECIES = {'homo_sapiens': 'hsa', 'mus_musculus': 'mmu', 'rattus_norvegicus': 'rno'}
MIRTRACE_PROTOCOL = config['filter']['trim']['mirtrace']['protocol']

def mirtrace_species():
    species = MIRTRACE_SPECIES.get(config['organism'])
    if species is None:
        raise ValueError('miRTrace QC supports {}, got organism {!r}'.format(', '.join(MIRTRACE_SPECIES), config['organism']))
    return species

rule mirtrace_config:
    input:
        expand(rules.merged_fastq_R1.output, sample=SAMPLES)
    output:
        join(MIRTRACE_DIR, 'mirtrace.config')
    params:
        adapter = config['adapter']
    run:
        with open(output[0], 'w') as fh:
            for sample, fastq in zip(SAMPLES, input):
                fh.write(','.join([os.path.abspath(fastq), sample, params.adapter]) + '\n')

rule qc_mirtrace:
    input:
        config = rules.mirtrace_config.output,
        fastq = expand(rules.merged_fastq_R1.output, sample=SAMPLES)
    output:
        results = join(MIRTRACE_DIR, 'mirtrace-results.json'),
        report = join(MIRTRACE_DIR, 'mirtrace-report.html'),
        length = join(MIRTRACE_DIR, 'mirtrace-stats-length.tsv'),
        contamination = join(MIRTRACE_DIR, 'mirtrace-stats-contamination_basic.tsv'),
        complexity = join(MIRTRACE_DIR, 'mirtrace-stats-mirna-complexity.tsv')
    params:
        outdir = MIRTRACE_DIR,
        protocol = MIRTRACE_PROTOCOL,
        species = mirtrace_species(),
        title = config.get('project_id', 'GCF-0000-00')
    threads:
        8
    container:
        'docker://' + config['docker']['mirtrace']
    shell:
        'mirtrace qc --force '
        '-s {params.species} -p {params.protocol} -c {input.config} '
        '-o {params.outdir} --title {params.title} -t {threads}'
