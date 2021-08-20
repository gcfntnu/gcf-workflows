#-*- mode:snakemake -*-
"""Indexes from fasta
"""
        
rule star_genome_index_gtf:
    input: 
        genome = join('{ref_dir}', 'fasta', '{prefix}.fa'),
        gtf = join('{ref_dir}', 'anno', 'genes.gtf')
    output:
        index = join('{ref_dir}', 'index', '{prefix}', 'star', 'r_{sjdbOverhang}', 'SA')
    params:
        index_dir =  join('{ref_dir}', 'index', '{prefix}', 'star', 'r_{sjdbOverhang}'),
        sjdbOverhang = '{sjdbOverhang}'
    threads:
        48
    singularity:
        'docker://' + config['docker']['star']
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.index_dir} '
        '--genomeFastaFiles {input.genome} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang {params.sjdbOverhang} '

rule star_genome_index:
    input: 
        genome = join('{ref_dir}', 'fasta', '{prefix}.fa')
    output:
        index = join('{ref_dir}', 'index', '{prefix}', 'star', 'SA')
    params:
        index_dir =  join('{ref_dir}', 'index', '{prefix}', 'star'),
        sjdbOverhang = '{sjdbOverhang}'
    threads:
        48
    singularity:
        'docker://' + config['docker']['star']
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.index_dir} '
        '--genomeFastaFiles {input.genome} '
        
rule hisat2_genome_index:
    input:
        genome = join('{ref_dir}', 'fasta', '{prefix}.fa')
    params:
        index_dir =  join('{ref_dir}', 'index', 'hisat2'),
        prefix = join('{ref_dir}', 'index', '{prefix}', 'hisat2', '{prefix}')
    output:
         join('{ref_dir}', 'index', '{prefix}', 'hisat2', '{prefix}'+'.1.ht2')
    threads:
        48
    log:
        join('{ref_dir}', 'logs', 'HISAT2.{prefix}.index.log')
    singularity:
        'docker://' + config['docker']['hisat2']
    shell:
        'hisat2-build '
        '-f '
        '-p {threads} '
        '{input.genome} '
        '{params.prefix} '

rule hisat2_splicesites:
    input:
        gtf = join('{ref_dir}', 'anno', 'genes.gtf')
    output:
        join('{ref_dir}', 'anno', 'splicesites.txt')
    singularity:
       'docker://' + config['docker']['hisat2'] 
    shell:
        'hisat2_extract_splice_sites.py {input.gtf} > {output}'
        
rule bwa_index:
    input:
        join('{ref_dir}', 'fasta', '{prefix}.fa')
    output:
        join('{ref_dir}', 'index', '{prefix}', 'bwa', '{prefix}.amb'),
        join('{ref_dir}', 'index', '{prefix}', 'bwa', '{prefix}.ann'),
        join('{ref_dir}', 'index', '{prefix}', 'bwa', '{prefix}.bwt'),
        join('{ref_dir}', 'index', '{prefix}', 'bwa', '{prefix}.pac'),
        join('{ref_dir}', 'index', '{prefix}', 'bwa', '{prefix}.sa')
    params:
        output_dir = join('{ref_dir}', 'index', '{prefix}', 'bwa'),
        input_dir = join('{ref_dir}', 'fasta')
    singularity:
        'docker://' + config['docker']['bwa_samtools']
    shell:
        """
        bwa index {input}
        mkdir -p {params.output_dir}
        mv {params.input_dir}/{wildcards.prefix}.fa.amb {params.output_dir}/{wildcards.prefix}.amb
        mv {params.input_dir}/{wildcards.prefix}.fa.ann {params.output_dir}/{wildcards.prefix}.ann
        mv {params.input_dir}/{wildcards.prefix}.fa.bwt {params.output_dir}/{wildcards.prefix}.bwt
        mv {params.input_dir}/{wildcards.prefix}.fa.pac {params.output_dir}/{wildcards.prefix}.pac
        mv {params.input_dir}/{wildcards.prefix}.fa.sa {params.output_dir}/{wildcards.prefix}.sa
        """

rule bowtie_index:
    input:
        join('{ref_dir}', 'fasta', '{prefix}.fa')
    output:
        join('{ref_dir}', 'index', '{prefix}', 'bowtie', '{prefix}.1.ebwt'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie', '{prefix}.2.ebwt'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie', '{prefix}.3.ebwt'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie', '{prefix}.4.ebwt'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie', '{prefix}.rev.1.ebwt'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie', '{prefix}.rev.2.ebwt')
    params:
        index = join('{ref_dir}', 'index', '{prefix}', 'bowtie', '{prefix}')
    singularity:
        'docker://' + config['docker']['bowtie']
    shell:
        'bowtie-build {input} {params.index}'

rule bowtie2_index:
    input:
        join('{ref_dir}', 'fasta', '{prefix}.fa')
    output:
        join('{ref_dir}', 'index', '{prefix}', 'bowtie2', '{prefix}.1.bt2'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie2', '{prefix}.2.bt2'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie2', '{prefix}.3.bt2'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie2', '{prefix}.4.bt2'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie2', '{prefix}.rev.1.bt2'),
        join('{ref_dir}', 'index', '{prefix}', 'bowtie2', '{prefix}.rev.2.bt2')
    params:
        index = join('{ref_dir}', 'index', '{prefix}', 'bowtie2', '{prefix}')
    threads:
        48
    singularity:
        'docker://' + config['docker']['bowtie2']
    shell:
        'bowtie2-build --threads {threads} {input} {params.index}'

rule bbmap_index:
    input:
         join('{ref_dir}', 'fasta', '{prefix}.fa')
    output:
        join('{ref_dir}', 'index', '{prefix}', 'bbmap', 'ref', '{prefix}', '1', 'info.txt')
    params:
        index = join('{ref_dir}', 'index', '{prefix}', 'bbmap')
    singularity:
        'docker://' + config['docker']['bbmap']
    shell:
        'bbmap.sh ref={input} path={params.index}'


rule salmon_selective_decoy:
    input:
        join('{ref_dir}', 'fasta', 'genome.fa')
    output:
        temp(join('{ref_dir}', 'decoys.txt'))
    shell:
        """
        grep "^>" {input} | cut -d " " -f1 > {output}
        sed -i -e 's/>//g' {output}
        """

rule salmon_selective_fasta:
    input:
        join('{ref_dir}', 'fasta', 'transcriptome.fa'),
        join('{ref_dir}', 'fasta', 'genome.fa')
    output:
        temp(join('{ref_dir}', 'fasta', 'gentrome.fa'))
    shell:
        'cat {input} > {output}'
        
ruleorder: salmon_index_selective > salmon_index


rule salmon_index_selective:
    input:
        fasta = rules.salmon_selective_fasta.output,
        decoy = join('{ref_dir}', 'decoys.txt')
    output:
        join('{ref_dir}', 'index', 'gentrome', 'salmon', 'info.json')
    params:
        out = join('{ref_dir}', 'index', 'gentrome', 'salmon')     
    threads:
        16
    priority:
        1
    singularity:
        'docker://' + config['docker']['salmon']
    shell:
        'salmon index '
        '-d {input.decoy} '
        '--threads {threads} '
        '--index {params.out} '
        '--transcripts {input.fasta}'
        
rule salmon_index:
    input:
        join('{ref_dir}', 'fasta', '{prefix}.fa')
    output:
        join('{ref_dir}', 'index', '{prefix}', 'salmon', 'info.json')
    params:
        out = join('{ref_dir}', 'index', '{prefix}', 'salmon')
    threads:
        16
    singularity:
        'docker://' + config['docker']['salmon']
    shell:
        'salmon index '
        '--threads {threads} '
        '--index {params.out} '
        '--transcripts {input}'


rule salmon_index_tximeta:
    input:
        fasta = join('{ref_dir}', 'fasta', '{prefix}.fa'),
        json = join('{ref_dir}', 'index', '{prefix}', 'salmon', 'info.json'),
        gtf = join('{ref_dir}', 'anno', 'genes.gtf')
    output:
        json = join('{ref_dir}', 'index', '{prefix}', 'salmon', 'tximeta.json')
    params:
        script = srcdir('scripts/make_linked_txome.R'),
        org = config['organism'],
        db = lambda x: 'GCF_' + config['db']['reference_db'],
        release = lambda x: DB_CONF['release'],
        assembly = lambda x: DB_CONF['assembly'],
        cache = join(EXT_CACHE, 'tximeta')
    threads:
        16
    singularity:
        'docker://' + config['docker']['tximport']
    shell:
        'Rscript {params.script} -v '
        '--index {input.json} '
        '--cachedir {params.cache} '
        '--transcriptome {input.fasta} '
        '--organism {params.org} '
        '--gtf {input.gtf} '
        '--source {params.db} '
        '--release {params.release} '
        '--assembly {params.assembly} '
        '-o {output.json} '
        '--verbose '
        
