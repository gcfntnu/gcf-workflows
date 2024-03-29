rule yanagi_preprocess:
    input:
        gtf = join(REF_DIR, 'anno', 'genes.gtf'),
        genome = join(REF_DIR, 'fasta', 'genome.fa')
    output:
        ex = join(QUANT_INTERIM, 'yanagi', 'preprocess', 'disjoint_bins.tsv'),
        tx = join(QUANT_INTERIM, 'yanagi', 'preprocess', 'txs2bins.tsv')
    params:
        outdir = join(QUANT_INTERIM, 'yanagi', 'preprocess')
    container:
        'docker://gcfntnu/yanagi:0.1'
    shell:
        """
        cd /opt/yanagi
        python yanagi.py preprocess -gtf {input.gtf} -fa {input.genome} -o {params.outdir}
        """

