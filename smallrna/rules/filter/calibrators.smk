#-*-mode:snakemake-*-

SPIKE_REF = config['filter']['spikein'].get('ref', '')

def get_spikein_index():
    ref = SPIKE_REF.lower()
    method = config['filter']['spikein'].get('quantifier')
    if not method in ['bowtie', 'bowtie2']:
        method = 'bowtie'
    if ref == 'ercc':
        return join(EXT_DIR, 'ERCC', 'index', 'ERCC92', method, 'ERCC92')
    elif ref == 'smallrna_calibrators':
        return join(EXT_DIR, 'spikein', 'index', 'smallrna_calibrators', method, 'smallrna_calibrators')
    else:
        return join(EXT_DIR, 'merged_spikein', 'index', 'merged_spikein', method, 'merged_spikein')

def get_spikein_fasta():
    ref = SPIKE_REF.lower()
    if ref == 'ercc':
        return join(EXT_DIR, 'ERCC', 'fasta', 'ERCC92.fa')
    elif ref == 'smallrna_calibrators':
        return join(EXT_DIR, 'spikein', 'fasta', 'smallrna_calibrators.fa')
    else:
        return join(EXT_DIR, 'spikein', 'fasta', 'merged_spikein.fa')


rule spikein_filter_bowtie:
    input:
        unpack(get_filtered_fastq),
        index = get_spikein_index() + '.1.ebwt'
    output:
        R1 = temp(join(FILTER_INTERIM, 'spikein', 'bowtie', '{sample}_R1.fastq')),
        counts = join(FILTER_INTERIM, 'spikein', 'bowtie', '{sample}.counts')
    container:
        'docker://' + config['docker']['bowtie_samtools']
    params:
        args = '-n 0 -k 1 -l 18 -q --best --norc -S ',
        index = get_spikein_index()
    threads:
        4
    log:
        log = join(FILTER_INTERIM, 'spikein', 'bowtie', '{sample}.filter.log'),
        error = join(FILTER_INTERIM, 'spikein', 'bowtie', '{sample}.filter.error')
    shell:
        'bowtie {params.index} '
        '{input.R1} '
        '--un {output.R1} '
        '-p {threads} '
        '{params.args} '
        '2>> {log.log} '
        '| samtools view -q 5 -S  - | cut -f3 | sort | uniq -c  > {output.counts} '
        '2>> {log.error} '


def _filter_get_calibrator_clean(wildcards):
    cal_quant = config['filter'].get('spikein', {}).get('quantifier', 'skip')
    if cal_quant == 'skip':
        return get_filtered_fastq(wildcards)
    if cal_quant is None and config['quant'].get('method') == 'unitas':
        cal_quant = 'unitas'
    if cal_quant == 'unitas':
        return get_filtered_fastq(wildcards)
    elif cal_quant in ['bowtie']:
        return join(FILTER_INTERIM, 'spikein', cal_quant, '{}_R1.fastq'.format(wildcards.sample))
    else:
        raise ValueError('calibrator quantifier option not valid: {}'.format(cal_quant))
