#-*- mode:snakemake -*- 

rule cutadapt_demultiplex:
    input:
        unpack(get_filtered_fastq)
    output:
        R1 = temp(join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R1.fastq')),
        R2 = temp(join(FILTER_INTERIM, 'cutadapt_demultiplex', '{sample}_R2.fastq')),
        log = join(FILTER_INTERIM, 'cutadapt_demultiplex', 'log', '{sample}_qiaseq_demultiplex.log'),
    params:
        script = srcdir('scripts/demultiplex_16s_its.py'),
        outdir = join(FILTER_INTERIM, 'cutadapt_demultiplex'),
        logdir = join(FILTER_INTERIM, 'cutadapt_demultiplex', 'log'),
        libkit = config["libprepkit"] + (" PE" if len(config["read_geometry"]) > 1 else " SE"),
        libprepconf = srcdir("../../../libprep.config")
    threads: 
        2
    singularity: 
        'docker://' + config['docker']['cutadapt']
    shell:
        '{params.script} '
        '--in1 {input.R1} '
        '--in2 {input.R2} '
        '--libkit "{params.libkit}" '
        '--libprepconfig {params.libprepconf} '
        '--output-dir {params.outdir} '
        '--log-dir {params.logdir} '

rule cutadapt_all:
    input:
        expand(join(FILTER_INTERIM, 'cutadapt_demultiplex', 'log', '{sample}_qiaseq_demultiplex.log'), sample=SAMPLES)
