#-*- mode: snakemake -*-
"""
"""
include: 
    join(GCFDB_DIR, 'illumina.db')
include: 
    join(GCFDB_DIR, 'univec.db')

rule illumina_fasta:
    input:
        join(EXT_DIR, 'illumina', 'adaptors', 'fasta', 'adaptors.fa'),
        join(EXT_DIR, 'univec_illumina', 'fasta', 'univec_illumina.fa')
    output:
        join(EXT_DIR, 'illumina', 'adaptors_ext', 'fasta', 'adaptors_ext.fa')
    shell:
        'cat {input} > {output}'


def fastq_screen_indexes(*args, **kw):
   """ Return name and index for screening databases.
   """
   # default human
   release = config['db'].get('ensembl', {}).get('release', '103')
   INDEXES = {'Human': join(EXT_DIR, 'ensembl', 'release-{}'.format(release), 'homo_sapiens', 'GRCh38', 'index', 'genome', 'bowtie2', 'genome.1.bt2')}
   org = config.get('organism')
   if org not in ['homo_sapiens', 'N/A', 'n/a']:
      if config['db'].get('reference_db') == 'ensembl':
         release = config['db']['ensembl']['release']
         assembly = config['db']['ensembl']['assembly']
         index = join(EXT_DIR, 'ensembl', 'release-{}'.format(release), org, assembly, 'index', 'genome', 'bowtie2', 'genome.1.bt2')
         INDEXES[org] = index
      elif config['db'].get('reference_db') == 'custom_reference':
         assembly = config['db']['ensembl']['assembly']
         index = join(EXT_DIR, 'custom_reference', org, assembly, 'index', 'genome', 'bowtie2', 'genome.1.bt2')
         INDEXES[org] = index
   # contamination (univec minus phiX/Illumina sequences)
   INDEXES['Contamination'] = join(EXT_DIR, 'univec_subset', 'index', 'univec_subset', 'bowtie2', 'univec_subset.1.bt2')
   # phiX
   INDEXES['PhiX'] = join(EXT_DIR, 'illumina', 'phix', 'index', 'phix', 'bowtie2', 'phix.1.bt2')
   # illumina specific sequences
   INDEXES['Adaptors/Primers'] = join(EXT_DIR, 'illumina', 'adaptors_ext', 'index', 'adaptors_ext', 'bowtie2', 'adaptors_ext.1.bt2')
   return INDEXES

rule fastq_screen_config:
    input:
        unpack(fastq_screen_indexes)
    output:
        fn = temp(os.path.abspath('fastq_screen.config'))
    run:
        import os
        with open(output.fn, 'w') as fh:
           for db_name, index in input.items():
              index_dir = os.path.dirname(index)
              base = os.path.basename(index).split('.')[0]
              index_base = os.path.join(index_dir, base)
              fh.write('DATABASE {} {}\n'.format(db_name.title(), index_base))


rule fastq_screen:
   input:
       unpack(get_filtered_fastq),
       config = rules.fastq_screen_config.output 
   output:
       join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt')
   params:
       args = '-q --force',
       subset = 400000,
       outdir = join(QC_INTERIM, 'fastq_screen')
   threads:
       4
   singularity:
       'docker://' + config['docker']['fastq-screen']
   shadow:
       'minimal'
   shell:
       'fastq_screen '
       '--aligner bowtie2 '
       '--threads {threads} '
       '--conf {input.config} '
       '--outdir . '
       '{params.args} '
       '{input.R1} '
       ' && mv {wildcards.sample}*_screen.txt {output} '
