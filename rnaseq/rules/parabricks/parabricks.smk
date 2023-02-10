#-*- mode: snakemake -*-
"""

Snakemake rules for aligning rna-seq fastq files to genome using the
Nvidia Clara Parabricks pipeline.

Dependencies
------------
STAR, https://github.com/alexdobin/STAR
Nvidia Clara Parabricks,

"""

PB_STAR_INTERIM = join(ALIGN_INTERIM, 'pb_star')

def pb_star_input_params(wildcards):
    """Multiple fastq files per sample workaround.

    STAR uses comma separated input when defining a read with mutiple fastq files.
    A comma separated string is not a valid input file in snakemake.
    A work around in snakemake is to build the comma separated input string as a params.fastq
    
    """
    fastq = get_filtered_fastq(wildcards)
    R1 = fastq['R1']
    R2 = fastq.get('R2', '')
    if isinstance(R1, (list, tuple)):
        R1 = ','.join(sorted(set(R1)))
    if isinstance(R2, (list, tuple)):
        R2 = ','.join(sorted(set(R2)))
    input_string = ' '.join([R1, R2])
    return input_string

def pb_star_genome_dir():
    sjdb_rlen = int(max(config['read_geometry'])) - 1
    return join(REF_DIR, 'index', 'genome', 'pb_star', 'r_{}'.format(sjdb_rlen))
        
rule pb_rna_align:
    input:
        unpack(get_filtered_fastq),
        genome = join(pb_star_genome_dir(), 'SA'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf'),
        genome_fasta = join(REF_DIR, 'fasta', 'genome.fa'),
    params:
        fastq = pb_star_input_params,
        prefix = lambda wildcards, output: output.bam.split('Aligned.sortedByCoord.out.bam')[0],
        genome_dir = pb_star_genome_dir(),
        #args = config['align']['star']['args'],
        sjdb_rlen = int(max(config['read_geometry'])) - 1,
        out_dir = join(PB_STAR_INTERIM, '{sample}'),
    output:
        bam = join(PB_STAR_INTERIM, '{sample}', '{sample}.Aligned.sortedByCoord.out.bam'),
        log = join(PB_STAR_INTERIM, '{sample}', '{sample}.STAR.log')
    threads:
       48
    singularity:
        'docker://' + config['docker']['parabricks']
    benchmark:
        'benchmarks/{sample}/pb_rna_fq2bam.txt'
    shell: 
        'pbrun rna_fq2bam '
        '--num-threads {threads} '
        '--genome-lib-dir {params.genome_dir} '
        '--ref {input.genome_fasta} '
        '--in-fq {params.fastq} '
        '--out-bam {output.bam} '
        '--output-dir {params.out_dir} '
        '--sjdb-overhang {params.sjdb_rlen} '
        '--read-files-command zcat '
        '--max-bam-sort-memory 100000000000 '
        '--logfile {output.log} '

rule pb_star_all:
    input:
        all_bams = expand(rules.pb_rna_align.output.bam, sample=SAMPLES)
