#-*- mode: snakemake -*-
"""
"""
from functools import lru_cache
from hashlib import file_digest
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
PROJECT_SCREEN_FASTA = None
if config.get('ref_dir') and ORG not in ('homo_sapiens', 'n/a') and config['db']['reference_db'] in (
        'ensembl', '10xgenomics', 'custom_reference', 'refseq', 'gencode'):
    PROJECT_SCREEN_FASTA = join(config['ref_dir'], 'fasta', 'genome.fa')
PROJECT_SCREEN_ALIAS = 'project_bowtie2_index'
PROJECT_SCREEN_CACHE = join(QC_INTERIM, 'fastq_screen', 'indexes')


@lru_cache(maxsize=1)
def project_screen_index():
    # Resolve on demand: workflows that do not request FastQ Screen need no extra FASTA read.
    with open(PROJECT_SCREEN_FASTA, 'rb') as fasta:
        digest = file_digest(fasta, 'sha256').hexdigest()
        # bowtie2-build 2.3.4.3 chooses the binary by FASTA byte size (small_index_max_size).
        extension = 'bt2l' if os.path.getsize(fasta.fileno()) > 4 * 1024**3 - 200 else 'bt2'
    return join(PROJECT_SCREEN_CACHE, digest, 'genome'), extension


def project_screen_setup(wildcards):
    if not PROJECT_SCREEN_FASTA:
        return ''
    prefix, extension = project_screen_index()
    return f'mkdir {PROJECT_SCREEN_ALIAS} && ln -s -- {quote(abspath(prefix))}.*.{extension} {PROJECT_SCREEN_ALIAS}/ && '


if PROJECT_SCREEN_FASTA:
    rule project_screen_bowtie2_index:
        input:
            PROJECT_SCREEN_FASTA
        output:
            bowtie2_index_files(join(PROJECT_SCREEN_CACHE, '{fasta_sha256}', 'genome'), '{index_ext}')
        params:
            script = join(GCFDB_DIR, 'scripts', 'build_project_screen_index.sh'),
            outdir = join(PROJECT_SCREEN_CACHE, '{fasta_sha256}')
        threads:
            48
        wildcard_constraints:
            fasta_sha256 = '[0-9a-f]{64}',
            index_ext = 'bt2l?'
        container:
            'docker://' + config['docker']['bowtie2']
        shell:
            'bash {params.script:q} {input:q} {wildcards.fasta_sha256:q} {wildcards.index_ext:q} {threads} {params.outdir:q}'


def fastq_screen_indexes(*args, **kw):
   """ Return name and index for screening databases.
   """
   # default human
   INDEXES = {'Human': join(HUMAN_ENSEMBL_DIR, 'index', 'genome', 'bowtie2', 'genome.1.bt2')}
   org = config.get('organism')
   if org not in ['homo_sapiens', 'N/A', 'n/a']:
      if PROJECT_SCREEN_FASTA:
         INDEXES[org] = bowtie2_index_files(*project_screen_index())
      elif config['db'].get('reference_db') in ['ensembl', '10xgenomics']:
         index = join(ENSEMBL_DIR, 'index', 'genome', 'bowtie2', 'genome.1.bt2')
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
        if PROJECT_SCREEN_FASTA and name == ORG:
            prefix = join(PROJECT_SCREEN_ALIAS, 'genome')
        else:
            prefix = abspath(str(index).removesuffix('.1.bt2'))
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
       project_setup = project_screen_setup,
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
