#-*- mode: snakemake -*-
"""
Snakemake rules for aligning rna-seq fastq files to genome using the
HISAT2 aligner.
"""

HISAT_INTERIM = join(ALIGN_INTERIM, 'hisat2')

PE = len(config['read_geometry']) > 1
KIT_STRAND = config['read_orientation']
if PE:
    if KIT_STRAND == 'forward':
        HISAT_STRAND = 'FR'
    else :
        HISAT_STRAND = 'RF'
else:
    if KIT_STRAND == 'forward':
        HISAT_STRAND = 'F'
    else :
        HISAT_STRAND = 'R'

def hisat2_input(wildcards):
    fastq = get_filtered_fastq(wildcards)
    sample = config['samples'][wildcards.sample]
    PE = sample.get('paired_end') or len(config['read_geometry']) > 1
    R1 = fastq['R1']
    R2 = fastq.get('R2', [])
    if isinstance(R1, str):
        R1 = [R1]
    if isinstance(R2, str):
        R2 = [R2]
    if PE:
        input = '-1 ' + ','.join(R1)
        R2 = fastq.get('R2', [])
        input += ' -2 ' + ','.join(R2)
    else:
        input = '-U ' + ','.join(R1)
    return input

rule hisat2_align:
    input:
        unpack(get_filtered_fastq),
        index = join(REF_DIR, 'index', 'genome', 'hisat2', 'genome' + '.1.ht2'),
        splices = join(REF_DIR, 'anno', 'splicesites.txt'),
        db_log = join(REF_DIR, 'logs', 'ensembl.dna.log')
    params:
        input = hisat2_input,
        index = lambda wildcards,input: input.index.split('.1.ht2')[0],
        strand = HISAT_STRAND,
        unaligned = join(HISAT_INTERIM, '{sample}.unaligned.fastq.gz')
    output:
        sam = temp(join(HISAT_INTERIM, '{sample}.out.sam'))
    threads:
        24
    container:
        'docker://' + config['docker']['hisat2']
    log:
        join('logs', '{sample}', '{sample}.hisat2.log')
    shell:
        'hisat2 '
        '--known-splicesite-infile {input.splices} '
        '--rna-strandness {params.strand} '
        '--summary-file {log} '
        '--new-summary '
        '-p {threads} '
        '--un-conc-gz {params.unaligned} '
        '-x {params.index} '
        '{params.input} '
        '-S {output.sam} '

rule hisat2_sorted_bam:
    input:
        sam = rules.hisat2_align.output.sam
    output:
        bam = temp(join(ALIGN_INTERIM, 'hisat2', '{sample}.Aligned.sortedByCoord.out.bam'))
    threads:
        4
    shell:
         'cat {input.sam} | sambamba view -S -f bam /dev/stdin | sambamba sort -m 24G -t {threads} -o {output.bam} /dev/stdin '

rule hisat2_all:
    input:
        expand(rules.hisat2_sorted_bam.output.bam, sample=SAMPLES)
    output:
        temp(touch('.hisat2.align.finalized'))
