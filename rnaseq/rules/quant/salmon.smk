#-*- mode: snakemake -*-
"""Snakemake rules for alignment free transcript and gene expression modelling using salmon.


Dependencies
------------
SALMON, https://combine-lab.github.io/salmon/

"""

SALMON_INTERIM = join(QUANT_INTERIM, 'salmon')
if not 'SALMON_INDEX_TYPE' in locals():
   SALMON_INDEX_TYPE = config.get('quant', {}).get('salmon', {}).get('index', 'transcriptome')
   
include:
    'tximport.smk'
    
if not 'REF_DIR' in locals():
   raise ValueError

   
if config.get('read_orientation') is None:
   LIB_TYPE = 'A'
else:
   PE = len(config['read_geometry']) > 1
   stranded = 'U' if config['read_orientation'] == 'unstranded' else 'S'
   if config['read_orientation'] == 'reverse':
      orientation = 'R'
   elif config['read_orientation'] == 'forward':
      orientation = 'F'
   else:
      orientation = ''
   LIB_TYPE = 'I' + stranded + orientation if PE else stranded + orientation


def get_fragment_length(wildcards):
    val = config['samples'][wildcards.sample].get('Fragment_Length', '180')
    if not val or val == 'NA':
        val = '180'
    return val
    
def get_fragment_sd(wildcards):
    val = config['samples'][wildcards.sample].get('Fragment_SD', '20')
    if not val or val == 'NA':
        val = '20'
    if int(val) < 20:
        val = '20'
    return val
    
def sample_args(wildcards):
    """Sample specific args
    """
    args = ''
    sample = config['samples'][wildcards.sample]
    PE = sample.get('paired_end') or len(config['read_geometry']) > 1
    fastq = get_filtered_fastq(wildcards)
    if isinstance(fastq['R1'], str):
        salmon_r1 = fastq['R1']
    else:
        salmon_r1 = ' '.join(fastq['R1'])
    if PE:
        if isinstance(fastq['R2'], str):
            salmon_r2 = fastq['R2']
        else:
            salmon_r2 = ' '.join(fastq['R2'])
        args += ' -1 {} -2 {} '.format(salmon_r1, salmon_r2)
        args += '--gcBias '
    else:
        args += ' -r {} '.format(salmon_r1)
        args += '--fldMean {} '.format(get_fragment_length(wildcards))
        args += '--fldSD {} '.format(get_fragment_sd(wildcards))
    return args


rule salmon_map:
    input:
        unpack(get_filtered_fastq),
        index = join(REF_DIR, 'index', SALMON_INDEX_TYPE, 'salmon', 'info.json'),
        gtf = join(REF_DIR, 'anno', 'genes.gtf')
    params:
        output = join(SALMON_INTERIM, '{sample}'),
        nboot = 100,
        lib_type = LIB_TYPE,
        index = join(REF_DIR, 'index', SALMON_INDEX_TYPE, 'salmon'),
        sample_args = sample_args,
        unmapped_args = '--writeUnmappedNames --writeMappings=unmapped.sam '
    threads:
        24
    container:
        'docker://' + config['docker']['salmon']
    output:
        quant = join(SALMON_INTERIM, '{sample}', 'quant.sf'),
        meta_log = join(SALMON_INTERIM, '{sample}', 'aux_info', 'meta_info.json'),
        dist_log = join(SALMON_INTERIM, '{sample}', 'libParams', 'flenDist.txt')
    shell:
        'salmon quant '
        '-i {params.index} '
        '-l {params.lib_type} '
        '-g {input.gtf} '
        '-p {threads} '
        '--numGibbsSamples {params.nboot} -d '
        '--validateMappings '
        '--rangeFactorizationBins 4 '
        '--seqBias '
        '-o {params.output} '
        '{params.sample_args} '
        
def unmapped_sample_args(wildcards):
    """Sample specific args for unmapped reads.
    """
    args = ''
    sample = config['samples'][wildcards.sample]
    PE = sample.get('paired_end') or len(config['read_geometry']) > 1
    fastq = get_filtered_fastq(wildcards)
    salmon_r1 = ' '.join(fastq['R1'])
    if PE:
        salmon_r2 = ' '.join(fastq['R2'])
        args += ' -1 {} -2 {} '.format(salmon_r1, salmon_r2)
    else:
        args += ' -r {} '.format(salmon_r1)
    return args

rule salmon_unmapped_fastq:
    input:
        unpack(get_filtered_fastq),
        unmapped = join(SALMON_INTERIM, '{sample}', 'unmapped.sam')
    params:
        args = unmapped_sample_args,
        prefix = join(SALMON_INTERIM, '{sample}', 'unmapped')
    container:
        'docker://' + config['docker']['salmon']
    output:
        join(SALMON_INTERIM, '{sample}', 'unmapped_1.fa.gz')
    shell:
        'salmontools extract-unmapped '
        '{params.args} '
        '-o {params.prefix} '


rule salmon_tximport:
    input:
        files = expand(join(SALMON_INTERIM, '{sample}', 'quant.sf'), sample=SAMPLES),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = srcdir('scripts/tximport.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        48
    output:
        rds = join(QUANT_INTERIM, 'salmon', 'tximport', '{}_salmon.rds'.format(SALMON_INDEX_TYPE))
    shell:
        'Rscript {params.script} {input.files} '
        '--txinfo {input.txinfo} '
        '-t salmon '
        '-o {output} '


rule salmon_tximeta:
    input:
        files = expand(join(SALMON_INTERIM, '{sample}', 'quant.sf'), sample=SAMPLES),
        sample_info = join(INTERIM_DIR, 'sample_info.tsv'),
        index = join(REF_DIR, 'index', '{prefix}', 'salmon', 'tximeta.json')
    params:
        script = srcdir('scripts/tximeta.R'),
        index_dir = join(REF_DIR, 'index', '{prefix}', 'salmon'),
        cache = join(EXT_CACHE, 'tximeta')
    threads:
        48
    container:
        'docker://' + config['docker']['tximport']
    output:
        rds = join(QUANT_INTERIM, 'salmon', 'tximeta', '{prefix}_salmon.rds')
    shell:
        'Rscript {params.script} {input.files} '
        '--sample-info {input.sample_info} '
        '--txome {input.index} '
        '--cachedir {params.cache} '
        '-t salmon '
        '--verbose '
        '-o {output} '

   
rule salmon_gene_anndata:
    input:
        gene_counts = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_counts.tsv'),
        gene_tpm = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_tpm.tsv'),
        gene_abundance = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_tpm_length_scaled.tsv'),
        gene_vst = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_vst.tsv'),
        gene_length = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_lengths.tsv')
    output:
        join(QUANT_INTERIM, 'salmon', 'gene_anndata.adh5')
    params:
        script = srcdir('scripts/create_anndata.py')
    threads:
        48
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--gene-counts {input.gene_counts} '
        '--tpm {input.gene_tpm} '
        '--vsn {input.gene_vst} '
        '--abundance {input.gene_abundance} '
        '--output {output} '

        
rule salmon_quant:
    input:
        gene_counts = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_counts.tsv'),
        gene_tpm = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_tpm.tsv'),
        gene_abundance = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_tpm_length_scaled.tsv'),
        gene_vst = join(QUANT_INTERIM, 'salmon', 'tximport', 'gene_vst.tsv'),
        tx_counts = join(QUANT_INTERIM, 'salmon', 'tximport', 'transcript_counts.tsv'),
        tx_tpm = join(QUANT_INTERIM, 'salmon', 'tximport', 'transcript_tpm.tsv')


rule terminus_group:
    input:
        join(SALMON_INTERIM, '{sample}', 'quant.sf')
    params:
        salmon_dir = join(SALMON_INTERIM, '{sample}'),
        outdir = join(QUANT_INTERIM, 'terminus'),
        m = 0.05,
        tolerance = 0.001
    container:
        'docker://' + config['docker']['salmon']
    threads:
        16       
    output:
        join(QUANT_INTERIM, 'terminus', '{sample}', 'groups.txt')
    shell:
        'terminus group -m {params.m} -d {params.salmon_dir} -o {params.outdir} '

rule terminus_collapse:
    input:
       salmon = expand(rules.terminus_group.input, sample=SAMPLES),
       groups = expand(rules.terminus_group.output, sample=SAMPLES)
    container:
        'docker://' + config['docker']['salmon']
    params:
        input_dirs = lambda wildcards, input: ' '.join([os.path.dirname(n) for n in input.salmon]),
        outdir = join(QUANT_INTERIM, 'terminus'),
        consensus_thresh = 0.5
    output:
        expand(join(QUANT_INTERIM, 'terminus', '{sample}', 'quant.sf'), sample=SAMPLES),
        expand(join(QUANT_INTERIM, 'terminus', '{sample}', 'clusters.txt'), sample=SAMPLES)
    threads:
        16  
    shell:
        'terminus collapse -c {params.consensus_thresh} -d {params.input_dirs} -o {params.outdir} '


rule terminus_tx2terminus:
    input:
        quant = join(SALMON_INTERIM, SAMPLES[0], 'quant.sf'),
        clusters = join(QUANT_INTERIM, 'terminus', SAMPLES[0], 'clusters.txt'),
        tx2gene = join(REF_DIR, 'anno', 'tx2gene.tsv')
    output:
        join(QUANT_INTERIM, 'terminus', 'tx2terminus.tsv')
    params:
        script = srcdir('scripts/extract_txp_group.py')
    threads:
        16
    shell:
        'python {params.script} {input.quant} {input.clusters} {input.tx2gene} {output}'
       
rule terminus_tximport:
    input:
        files = expand(join(SALMON_INTERIM, '{sample}', 'quant.sf'), sample=SAMPLES),
        txinfo = join(QUANT_INTERIM, 'terminus', 'tx2terminus.tsv')
    params:
        script = srcdir('scripts/tximport.R')
    container:
        'docker://' + config['docker']['tximport']
    output:
        rds = join(QUANT_INTERIM, 'terminus', 'tximport', '{}_terminus.rds'.format(SALMON_INDEX_TYPE))
    threads:
        48
    shell:
        'Rscript {params.script} {input.files} '
        '--txinfo {input.txinfo} '
        '-t terminus '
        '-o {output} '

rule terminus_info:
    input:
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv'),
        tx2terminus = join(QUANT_INTERIM, 'terminus', 'tx2terminus.tsv')
    params:
        script = srcdir('scripts/terminus_anno.py')
    output:
        join(QUANT_INTERIM, 'terminus', 'terminus_info.tsv')
    threads:
        16
    shell:
        'python {params.script} {input} {output}'
