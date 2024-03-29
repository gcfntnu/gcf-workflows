#-*- mode:snakemake -*-
"""QIMME2 classifiers


"""
from os.path import join
# libprep

libprep_conf_fn = src_gcf('../libprep.config')
if os.path.exists(libprep_conf_fn):
    with open(libprep_conf_fn) as fh:
        ALL_LIBPREP  = yaml.load(fh, Loader=Loader) or {}
LIBPREPS_16S = ['16S Metagenomic Sequencing Library Prep PE', 'QIAseq 16S ITS Region Panels PE']

DB = ['silva']
GG_RELEASES = ['13_8']
SILVA_RELEASES = ['132']
UNITE_RELEASES = ['8_2']
LEVELS = ['99']
UNITE_LEVELS = ['dynamic']

QIASEQ_16S = []
QIASEQ_ITS = []
for k, v in ALL_LIBPREP['QIAseq 16S ITS Region Panels PE']['db']['primers'].items():
    if k.startswith('V'):
        QIASEQ_16S.append(v)
    else:
        QIASEQ_ITS.append(v)
        
QIASEQ_16S.append('universal')
#ILMN_16S = ['V4']
ILMN_16S = list(ALL_LIBPREP['16S Metagenomic Sequencing Library Prep PE']['db']['primers'].values())
ILMN_ITS = list(ALL_LIBPREP['ITS Low Input GCF Custom PE']['db']['primers'].values())

ruleorder: prebuild_classifiers > build_q2_classifer

def get_primers(wildcards):
    for k, v in ALL_LIBPREP.items():
        if v['name'] == wildcards.libprep or k == wildcards.libprep.replace(' ', '_'):
            p = v['primers'][wildcards.region]
            return [p['forward'], p['reverse']]
    raise ValueError

rule prebuild_classifiers:
    params:
        url_1 = 'https://data.qiime2.org/2019.10/common/silva-132-99-nb-classifier.qza',
        url_2 = 'https://data.qiime2.org/2019.10/common/gg-13-8-99-nb-classifier.qza',
        proxy = config.get('proxy', {}).get('wget', ''),
    output:
       out_1 = join(EXT_DIR, 'silva', '132', 'qiime2', 'classifiers', '99_universal.qza'),
       out_2 = join(EXT_DIR, 'greengenes', '13_8', 'qiime2', 'classifiers', '99_universal.qza')
    shell:
        """
        wget {params.proxy} {params.url_1} -O- > {output.out_1}
        wget {params.proxy} {params.url_2} -O- > {output.out_2}
        """

rule build_q2_classifer:
    input:
        ref = join(EXT_DIR, '{db}', '{release}', 'fasta', '{level}', 'otus.fa'),
        taxa = join(EXT_DIR, '{db}', '{release}', 'anno', '{level}', 'otus_taxonomy.txt')
    params:
        script = src_gcf('scripts/train_classifier.py')
    output:
        join(EXT_DIR, '{db}', '{release}', 'qiime2', 'classifiers', '{level}_{fwd}-{rev}.qza')
    container:
        'docker://' + config['docker']['qiime2']
    threads:
        24
    shell:
        'python {params.script} '
        '--fasta {input.ref} '
        '--taxa {input.taxa} '
        '--f-primer {wildcards.fwd} '
        '--r-primer {wildcards.rev} '
        '--min-len 30 '
        '--output {output} '
        '--threads {threads} '

rule export_classifier:
    input:
        join(EXT_DIR, '{db}', '{release}', 'qiime2', 'classifiers', '{level}_{region}.qza')
    output:
        join(EXT_DIR, '{db}', '{release}', 'qiime2', 'classifiers', 'export', '{level}_{region}', 'sklearn_pipeline.tar')
    params:
        outdir = join(EXT_DIR, '{db}', '{release}', 'qiime2', 'classifiers', 'export', '{level}_{region}')
    container:
        'docker://' + config['docker']['qiime2']
    shell:
        'qiime tools export --input-path {input} --output-path {params.outdir}'
        
rule q2_classifiers_all:
    input:
        expand(join(EXT_DIR, 'silva', '{release}', 'qiime2', 'classifiers', 'export', '{level}_{region}', 'sklearn_pipeline.tar'),
               release=SILVA_RELEASES, region=QIASEQ_16S, level=LEVELS),
        expand(join(EXT_DIR, 'greengenes', '{release}', 'qiime2', 'classifiers', 'export', '{level}_{region}', 'sklearn_pipeline.tar'),
               release=GG_RELEASES, region=QIASEQ_16S, level=LEVELS),
        expand(join(EXT_DIR, 'unite', '{release}', 'qiime2', 'classifiers', 'export', '{level}_{region}', 'sklearn_pipeline.tar'),
               release=UNITE_RELEASES, region=QIASEQ_ITS, level=UNITE_LEVELS),
        expand(join(EXT_DIR, 'unite', '{release}', 'qiime2', 'classifiers', 'export', '{level}_{region}', 'sklearn_pipeline.tar'),
                release=UNITE_RELEASES, region=ILMN_ITS, level=UNITE_LEVELS)

