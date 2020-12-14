#-*- mode:snakemake -*-


rule seqbuster_tabular:
    input:
        join(FILTER_INTERIM, 'mirtrace', 'qc_passed_reads.all.collapsed', '{sample}_R1.fasta')
    output:
        join(FILTER_INTERIM, 'seqbuster', '{sample}_R1.tsv')
    params:
        script = srcdir('scripts/seqbuster_tabular.py')
    shell:
        'python {params.script} {input} > {output}'

rule test_tab_all:
    input:
        expand(join(FILTER_INTERIM, 'seqbuster', '{sample}_R1.tsv'), sample=SAMPLES)


rule seqbuster_align:
    input:
        fasta = join(FILTER_INTERIM, 'seqbuster', '{sample}_R1.tsv'),
        db = join(MIRBASE_EXT, 'seqbuster')
    params:
        species = MIR_ORG[config['organism']],
        prefix = join(QUANT_INTERIM, 'seqbuster', '{sample}'),
        args = '-sub 1 -trim 3 -add 3'
    output:
        join(QUANT_INTERIM, 'seqbuster', '{sample}.mirna')
    singularity:
        'docker://' + config['docker']['seqbuster']
    threads:
        2
    shell:
        'miraligner '
        '{params.args} '
        '-s {params.species} '
        '-i {input.fasta} '
        '-db {input.db}  '
        '-o {params.prefix} '

rule seqbuster_mirtop:
    input:
       aligned = join(QUANT_INTERIM, 'seqbuster', '{sample}.mirna'),
       hairpin = join(MIRBASE_EXT, 'seqbuster', 'hairpin.fa'),
       gtf = join(MIRBASE_EXT, 'seqbuster', 'genes.gff')
    params:
        outdir = join(QUANT_INTERIM, 'seqbuster'),
        species = MIR_ORG[config['organism']]
    output:
        join(QUANT_INTERIM, 'seqbuster', '{sample}.gff')
    singularity:
        'docker://' + config['docker']['seqbuster']
    threads:
        2
    shell:
        'mirtop gff '
        '--format seqbuster '
        '--sps {params.species} '
        '--hairpin {input.hairpin}'
        '--gtf {input.gtf} '
        '-o {params.outdir} '
        '{input} '

rule seqbuster_isomir:
    input:
        aligned = expand(join(QUANT_INTERIM, 'seqbuster', '{sample}.mirna'), sample=SAMPLES),
        sample_info = 'data/tmp/sample_info.tsv'
    params:
        script = srcdir('scripts/seqbuster_count.R'),
        args = '--output-expression --output-psi ',
        outdir = join(QUANT_INTERIM, 'seqbuster', 'isomirs')
    threads:
        48
    output:
       join(QUANT_INTERIM, 'seqbuster', 'isomirs', 'mir_counts.txt'),
       join(QUANT_INTERIM, 'seqbuster', 'isomirs', 'isomir_counts.txt'),
       join(QUANT_INTERIM, 'seqbuster', 'isomirs', 'mir_vst.txt'),
       join(QUANT_INTERIM, 'seqbuster', 'isomirs', 'isomir_vst.txt'),
       join(QUANT_INTERIM, 'seqbuster', 'isomirs', 'mir_anno.txt'),
       join(QUANT_INTERIM, 'seqbuster', 'isomirs', 'isomir_anno.txt')
    singularity:
        'docker://' + config['docker']['seqbuster']    
    shell:
        'Rscript {params.script} '
        '--output-dir {params.outdir} '
        '--sample-info {input.sample_info} '
        '{params.args} '
        '{input.aligned} '
        

        
rule seqbuster_all:
    input:
        
