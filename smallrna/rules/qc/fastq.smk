#-*- mode: snakemake -*-
"""Fastq qc rules on trimmed fastq files
"""

MIRTRACE_DIR = join(INTERIM_DIR, 'smallrna', 'qc', 'mirtrace')

rule mirtrace_config:
    input:
        expand(join(FILTER_INTERIM, 'fastq', '{sample}_R1.fastq'), sample=SAMPLES)
    output:
        temp('.filter_mirtrace.conf')
    params:
        adapter = config['adapter']
    run:
        with open(output[0], 'w') as fh:
            for k, v in config['samples'].items():
                fn = join(config['fastq_dir'], v['R1'])
                els = [fn, k, params.adapter]
                fh.write(','.join(els) + '\n')            

def mirtrace_species(*args):
    org_map = {'homo_sapiens': 'hsa',
               'mus_musculus': 'mmu',
               'rattus_norvegicus': 'rno'}
    org = org_map.get(config['organism'])
    if org is None:
        raise ValueError('mirtrace qc support onle for human, mouse or rat. Got: {}'.format(config['organism']))
    return org

rule qc_mirtrace:
    input:
        config = '.filter_mirtrace.conf',
        samples = expand(join(FILTER_INTERIM, 'fastq', '{sample}.fastq'), sample=SAMPLES)
    params:
        outdir = MIRTRACE_DIR,
        protocol = config['filter']['trim']['mirtrace']['protocol'],
        species = mirtrace_species(),
        title = config.get('project_id', 'GCF-0000-00')
    container:
        'docker://' + config['docker']['mirtrace']
    threads:
        8
    output:
        qc = join(MIRTRACE_DIR, 'mirtrace-results.json'),
        length = join(MIRTRACE_DIR, 'mirtrace-stats-length.tsv')
    shell:
        'mirtrace qc '
        '-s {params.species} '
        '-c {input.config} '
        '-o {params.outdir} '
        '-p {params.protocol} '
        '--title {params.title} '
        '-t {threads} '
