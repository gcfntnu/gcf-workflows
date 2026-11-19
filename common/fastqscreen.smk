#-*- mode: snakemake -*-
"""
"""
from shlex import quote

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


# FastQ Screen treats whitespace as a field separator in DATABASE records. Each screening
# job uses its own shadow directory, so its project index can have a fixed relative alias.
PROJECT_SCREEN_DB = None
if config.get('ref_dir') and ORG not in ('homo_sapiens', 'n/a') and config['db']['reference_db'] in (
        'ensembl', '10xgenomics', 'custom_reference', 'refseq', 'gencode'):
    PROJECT_SCREEN_DB = join(config['ref_dir'], 'index', 'genome', 'bowtie2')
PROJECT_SCREEN_ALIAS = 'project_bowtie2_index'
PROJECT_SCREEN_SETUP = f'ln -s -- {quote(abspath(PROJECT_SCREEN_DB))} {PROJECT_SCREEN_ALIAS} && ' if PROJECT_SCREEN_DB else ''


def fastq_screen_indexes(*args, **kw):
   """ Return name and index for screening databases.
   """
   # default human
   INDEXES = {'Human': join(HUMAN_ENSEMBL_DIR, 'index', 'genome', 'bowtie2', 'genome.1.bt2')}
   org = config.get('organism')
   if org not in ['homo_sapiens', 'N/A', 'n/a']:
      if PROJECT_SCREEN_DB:
         INDEXES[org] = join(PROJECT_SCREEN_DB, 'genome.1.bt2')
      elif config['db'].get('reference_db') in ['ensembl', '10xgenomics']:
         release = config['db']['ensembl']['release']
         assembly = config['db']['ensembl']['assembly']
         index = join(EXT_DIR, 'ensembl', 'release-{}'.format(release), org, assembly, 'index', 'genome', 'bowtie2', 'genome.1.bt2')
         INDEXES[org] = index
      elif config['db'].get('reference_db') == 'custom_reference':
         assembly = config['db']['custom_reference'].get('assembly', '')
         index = join(EXT_DIR, 'custom_reference', org, assembly, 'index', 'genome', 'bowtie2', 'genome.1.bt2')
         INDEXES[org] = index
   # contamination (univec minus phiX/Illumina sequences)
   INDEXES['Contamination'] = join(EXT_DIR, 'univec_subset', 'index', 'univec_subset', 'bowtie2', 'univec_subset.1.bt2')
   # phiX
   INDEXES['PhiX'] = join(EXT_DIR, 'illumina', 'phix', 'index', 'phix', 'bowtie2', 'phix.1.bt2')
   # illumina specific sequences
   INDEXES['Adaptors/Primers'] = join(EXT_DIR, 'illumina', 'adaptors_ext', 'index', 'adaptors_ext', 'bowtie2', 'adaptors_ext.1.bt2')
   return INDEXES

def fastq_screen_databases(wildcards, input):
    records = []
    for name, index in input.items():
        prefix = abspath(str(index).removesuffix('.1.bt2'))
        if PROJECT_SCREEN_DB and name == ORG:
            prefix = join(PROJECT_SCREEN_ALIAS, 'genome')
        records.append(f'DATABASE {name.title()} {prefix}')
    return records


rule fastq_screen_config:
    input:
        unpack(fastq_screen_indexes)
    output:
        fn = temp(os.path.abspath('fastq_screen.config'))
    params:
        databases = fastq_screen_databases
    shell:
        "printf '%s\\n' {params.databases:q} > {output.fn:q}"


rule fastq_screen:
   input:
       unpack(get_filtered_fastq),
       config = rules.fastq_screen_config.output 
   output:
       join(QC_INTERIM, 'fastq_screen', '{sample}_screen.txt')
   params:
       project_setup = PROJECT_SCREEN_SETUP,
       args = '-q --force',
       subset = 400000,
       outdir = join(QC_INTERIM, 'fastq_screen')
   threads:
       4
   container:
       'docker://' + config['docker']['fastq-screen']
   shadow:
       'minimal'
   shell:
       '{params.project_setup}fastq_screen '
       '--aligner bowtie2 '
       '--threads {threads} '
       '--conf {input.config} '
       '--outdir . '
       '{params.args} '
       '{input.R1} '
       ' && mv {wildcards.sample}*_screen.txt {output} '
