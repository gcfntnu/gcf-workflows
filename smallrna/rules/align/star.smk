#-*- mode: snakemake -*-
"""

Snakemake rules for aligning rna-seq fastq files to genome using the
STAR aligner.

This is a single pass alignment with known reference and gene model. The
output is coordinate sorted bam files with marked duplicates and the
companion index (.bai) file.


Dependencies
------------
STAR, https://github.com/alexdobin/STAR
SAMBAMBA,

"""

STAR_INTERIM = join(ALIGN_INTERIM, 'star')

def star_input_params(wildcards):
    """Multiple fastq files per sample workaround.

    STAR uses comma separated input when defining a read with mutiple fastq files.
    A comma separated string is not a valid input file in snakemake.
    A work around in snakemake is to build the comma separated input string as a params.fastq
    
    """
    R1 = expand(get_filtered_fastq(wildcards)['R1'], sample=wildcards.sample)
    if isinstance(R1, str):
        R1 = [R1]
    input_string = ' '.join(R1)
    return input_string

def star_genome_dir():
    sjdb_rlen = int(max(config['read_geometry'])) - 1
    return join(REF_DIR, 'index', 'genome', 'star', 'r_{}'.format(sjdb_rlen))
        
rule star_align:
    input:
        unpack(get_filtered_fastq),
        genome = join(star_genome_dir(), 'SA'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    params:
        fastq = star_input_params,
        prefix = lambda wildcards, output: output.bam.split('Aligned.sortedByCoord.out.bam')[0],
        genome_dir = star_genome_dir(),
        args = config['align']['star']['args']
    output:
        bam = temp(join(STAR_INTERIM, '{sample}.Aligned.sortedByCoord.out.bam')),
        log = join(STAR_INTERIM, '{sample}.Log.final.out')
    threads:
       16
    container:
        'docker://' + config['docker']['star']
    shell: 
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outFileNamePrefix {params.prefix} '
        '--readFilesIn {params.fastq} '
        '{params.args} '

rule star_all:
    input:
        all_bams = expand(rules.star_align.output.bam, sample=SAMPLES)
    params:
        genome_dir = star_genome_dir()
    output:
        temp(touch('.star.align.finalized'))
    container:
        'docker://' + config['docker']['star']
    priority:
        0
    shell:
        'STAR '
        '--genomeDir {params.genome_dir} '
        '--genomeLoad Remove '
        '--outFileNamePrefix /tmp/foo '
        '|| echo "NO STAR SHARED MEMORY" '
