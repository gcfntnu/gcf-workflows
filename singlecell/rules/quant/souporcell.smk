"""souporcell
souporcell: Robust clustering of single cell RNAseq by genotype and ambient RNA inference without reference genotypes

https://www.biorxiv.org/content/10.1101/699637v2.full.pdf
https://github.com/wheaton5/souporcell
"""

SPC_INTERIM = join(QUANT_INTERIM, 'souporcell')

rule spc_common_snps:
    output:
        join(SPC_INTERIM, 'spc_common_snvs_hg38')
    params:
        url = 'https://data.genomicsresearch.org/Projects/scSplit/CommonSNVs/common_snvs_hg38.tar.gz'
    shell:
        """
        wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_" -O {output} && rm -rf /tmp/cookies.txt
        """


rule souporcell_barcodes:
    input:
        rules.cellranger_quant.output.filt_barcodes
    output:
        temp('{sample}_barcodes.tsv')
    shell:
        'gunzip -c {input} > {output}'

rule souporcell_run:
    input:
        bam = rules.cellranger_quant.output.bam,
        barcodes = rules.souporcell_barcodes.output,
        genome = join(CR_REF_DIR, 'fasta', 'genome.fa')
    params:
        k = config.get('souporcell',{}).get('num_mixed', 4),
        out_dir = join(SPC_INTERIM, '{sample}')
    output:
        join(SPC_INTERIM, '{sample}', 'cluster_genotypes.vcf')
    container:
        'docker://gcfntnu/souporcell:2.0'
    threads:
        24
    shell:
        'souporcell_pipeline.py -'
        'i {input.bam} '
        '-b {input.barcodes} '
        '-f {input.genome} '
        '-t {threads} '
        '-o {params.out_dir} '
        '-k {params.k} '

rule souporcell_all:
    input:
        expand(rules.souporcell_run.output, sample=SAMPLES)
