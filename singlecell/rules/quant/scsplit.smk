"""scSplit
Genotype-free demultiplexing of pooled single-cell RNA-seq, using a
hidden state model for identifying genetically distinct samples within
a mixed population.

https://doi.org/10.1186/s13059-019-1852-7
https://github.com/jon-xu/scSplit
"""

SCSPLIT_INTERIM = join(QUANT_INTERIM, 'scsplit')


def optical_dup_args(*args, **kw):
    args = 'TAGGING_POLICY=All '
    machine = config.get('machine', 'NextSeq 500')
    if machine == 'NextSeq 500':
        args += 'OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 '
    elif machine == 'HiSeq 2500':
        args += 'OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 '
    else:
        args = ''
    if config.get('remove_optical_duplicates', True):
        args += 'REMOVE_SEQUENCING_DUPLICATES=TRUE '
    return args

rule picard_mark_duplicates:
    input:
        bam = join(STAR_INTERIM, '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        bam = join(SCSPLIT_INTERIM, '{sample}', 'mrkdup.bam'),
        metrics = join(SCSPLIT_INTERIM, 'logs','{sample}.rmdup.metrics')
    log:
        join(SCSPLIT_INTERIM, 'logs','{sample}.picard.rmdup.log')
    params:
        java_opt="-Xms4g -Xmx4g ",
        dup_args = optical_dup_args()
    threads:
        8
    container:
        'docker://' + config['docker']['picard_gatk']
    shell:
        'picard MarkDuplicates '
        '{params.java_opt} '
        'INPUT={input.bam} '
        'OUTPUT={output.bam} '
        'METRICS_FILE={output.metrics} '
        'VALIDATION_STRINGENCY=SILENT '
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
        '2> {log} '

rule clean_bam:
    input:
        bam = join(STAR_INTERIM, '{sample}', 'Aligned.sortedByCoord.out.bam'),
    output:
        bam = join(SCSPLIT_INTERIM, '{sample}', 'filtered.bam')
    container:
        'docker://' + config['docker']['bowtie2_samtools']
    params:
        '-q 10 -F 3844 '
    shell:
        """
        samtools view -S -b {params} {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule umitools_dedup:
    input:
        bam = join(SCSPLIT_INTERIM, '{sample}', 'filtered.bam')
    output:
        bam = join(SCSPLIT_INTERIM, '{sample}', 'dedup.bam')
    container:
        'docker://' + config['docker']['umi_tools']
    log:
        join(SCSPLIT_INTERIM, 'logs', '{sample}.umitools.dedup.txt')
    params:
        args = '--extract-umi-method=tag --cell-tag=CB --umi-tag=UB --method=unique '
    shell:
        'umi_tools dedup --stdin={input.bam} --output-stats={log}  -S {output} {params.args} '

rule sort_bam:
    input:
        bam = join(SCSPLIT_INTERIM, '{sample}', 'dedup.bam')
    output:
        bam = join(SCSPLIT_INTERIM, '{sample}', 'sorted.bam')
    container:
        'docker://' + config['docker']['bowtie2_samtools']
    shell:
        'samtools sort {input} -o {output}'

rule index_bam:
    input:
        bam = join(SCSPLIT_INTERIM, '{sample}', 'sorted.bam') 
    output:
        bai = join(SCSPLIT_INTERIM, '{sample}', 'sorted.bam.bai')
    container:
        'docker://' + config['docker']['bowtie2_samtools']
    shell:
        'samtools index {input}'

rule freebayes_varcall:
    input:
        genome = join(EXT_DIR, 'ensembl/homo_sapiens/release-100/GRCh38/fasta', 'genome.fa'),
        genome_index = join(EXT_DIR, 'ensembl/homo_sapiens/release-100/GRCh38/fasta', 'genome.fa.fai'),
        bam = join(SCSPLIT_INTERIM, '{sample}', 'sorted.bam'),
        bai = join(SCSPLIT_INTERIM, '{sample}', 'sorted.bam.bai')
    threads:
        24
    output:
        vcf = join(SCSPLIT_INTERIM, '{sample}', 'var.vcf')
    shell:
        'freebayes-parallel <(fasta_generate_regions.py {input.genome_index} 100000) 24 -f {input.genome} -iXu -C 2 -q 1 {input.bam} > {output} '

rule scsplit_call:
    input:
        expand(rules.freebayes_varcall.output, sample=SAMPLES)

rule common_snps:
    output:
        join(SCSPLIT_INTERIM, 'common_snvs_hg38')
    params:
        url = 'https://data.genomicsresearch.org/Projects/scSplit/CommonSNVs/common_snvs_hg38.tar.gz'
    shell:
        'wget {params.url} -qO- | tar xvz -C {SCSPLIT_INTERIM} '
        
rule scsplit_count:
    input:
       vcf = join(SCSPLIT_INTERIM, '{sample}', 'var.vcf'),
       bam = join(SCSPLIT_INTERIM, '{sample}', 'sorted.bam'),
       common_snvs = join(SCSPLIT_INTERIM, 'common_snvs_hg38'),
       whitelist = join(STAR_INTERIM, '{sample}', 'Solo.out', 'Gene', 'filtered', 'barcodes.tsv')
    output:
        ref = join(SCSPLIT_INTERIM, '{sample}', 'ref_filtered.csv'),
        alt = join(SCSPLIT_INTERIM, '{sample}', 'alt_filtered.csv')
    params:
        out_dir = join(SCSPLIT_INTERIM, '{sample}')
    threads:
        24
    container:
        'docker://gcfntnu/scsplit:1.0.8'
    shell:
        'scSplit count -v {input.vcf} -i {input.bam} -c {input.common_snvs} -b {input.whitelist} '
        '-r ref_filtered.csv -a alt_filtered.csv -o {params.out_dir} '
        
rule scsplit_run:
    input:
        ref = join(SCSPLIT_INTERIM, '{sample}', 'ref_filtered.csv'),
        alt = join(SCSPLIT_INTERIM, '{sample}', 'alt_filtered.csv')
    params:
        out_dir = join(SCSPLIT_INTERIM, '{sample}'),
        n = config.get('scsplit',{}).get('num_mixed', 4),
        doublet_rate = 0.1
    output:
        res = join(SCSPLIT_INTERIM, '{sample}', 'scSplit_result.csv'),
        var = join(SCSPLIT_INTERIM, '{sample}', 'scSplit_dist_variants.txt'),
        dist = join(SCSPLIT_INTERIM, '{sample}', 'scSplit_dist_matrix.csv'),
        PA = join(SCSPLIT_INTERIM, '{sample}', 'scSplit_PA_matrix.csv'),
        P_s_c = join(SCSPLIT_INTERIM, '{sample}', 'scSplit_P_s_c.csv'),
        log = join(SCSPLIT_INTERIM, '{sample}', 'scSplit.log')
    container:
        'docker://gcfntnu/scsplit:1.0.8'
    threads:
        24
    shell:
        'scSplit run -r {input.ref} -a {input.alt} -n {params.n} -o {params.out_dir} -d {params.doublet_rate} '
        
rule scsplit_genotype:
    input:
        ref = join(SCSPLIT_INTERIM, '{sample}', 'ref_filtered.csv'),
        alt = join(SCSPLIT_INTERIM, '{sample}', 'alt_filtered.csv'),
        P_s_c = join(SCSPLIT_INTERIM, '{sample}', 'scSplit_P_s_c.csv')
    params:
        out_dir = join(SCSPLIT_INTERIM, '{sample}')
    output:
        join(SCSPLIT_INTERIM, '{sample}', 'scSplit.vcf')
    container:
        'docker://gcfntnu/scsplit:1.0.8'
    shell:
        'scSplit genotype -r {input.ref} -a {input.alt} -o {params.out_dir} -p {input.P_s_c} '        


rule scsplit_scanpy:
    input:
        h5ad = join(QUANT_INTERIM, 'aggregate', 'star', 'scanpy', 'all_samples_aggr.h5ad'),
        res = expand(join(SCSPLIT_INTERIM, '{sample}', 'scSplit_result.csv'), sample=SAMPLES)
    output:
        h5ad = join(SCSPLIT_INTERIM, 'scannpy_scsplit.h5ad')
    params:
        script = workflow.source_path('scripts/scsplit_scanpy.py')
    shell:
        'python {params.script} '
        '--scsplit {input.res} '
        '--h5ad {input.h5ad} '
        
