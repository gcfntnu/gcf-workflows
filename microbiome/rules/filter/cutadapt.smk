#-*- mode:snakemake -*- 

rule cutadapt_demultiplex:
    input:
        unpack(get_raw_fastq)
    output:
        R1 = temp(join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R1.fastq')),
        R2 = temp(join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R2.fastq')),
        log = join(FILTER_INTERIM, 'cutadapt_demultiplex', 'log', '{sample}_qiaseq_demultiplex.log'),
    params:
        script = src_gcf('scripts/demultiplex_16s_its.py'),
        outdir = join(FILTER_INTERIM, 'cutadapt_demultiplex'),
        logdir = join(FILTER_INTERIM, 'cutadapt_demultiplex', 'log'),
        libkit = config["libprepkit"] + (" PE" if len(config["read_geometry"]) > 1 else " SE"),
        libprepconf = src_gcf("../../../libprep.config")
    threads: 
        2
    container: 
        'docker://' + config['docker']['cutadapt']
    shell:
        '{params.script:q} '
        '--in1 {input.R1:q} '
        '--in2 {input.R2:q} '
        '--libkit {params.libkit:q} '
        '--libprepconfig {params.libprepconf:q} '
        '--output-dir {params.outdir:q} '
        '--log-dir {params.logdir:q} '

rule cutadapt_all:
    input:
        expand(join(FILTER_INTERIM, 'cutadapt_demultiplex', 'log', '{sample}_qiaseq_demultiplex.log'), sample=SAMPLES)
