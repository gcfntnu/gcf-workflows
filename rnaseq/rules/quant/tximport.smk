if not 'SALMON_INDEX_TYPE' in locals():
   SALMON_INDEX_TYPE = config.get('quant', {}).get('salmon', {}).get('index', 'transciptome')

rule tximport_gene_counts:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_counts.tsv')
    threads:
        8
    shell:
        'Rscript {params.script} '
        '--type gene '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_gene_lengths:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_lengths.tsv')
    threads:
        8
    shell:
        'Rscript {params.script} '
        '--type gene_length '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_gene_tpm:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_tpm.tsv')
    shell:
        'Rscript {params.script} '
        '--type gene_tpm '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '

rule tximport_gene_tpm_scaled:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_tpm_scaled.tsv')
    shell:
        'Rscript {params.script} '
        '--type gene_tpm '
        '--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_gene_tpm_length_scaled:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_tpm_length_scaled.tsv')
    shell:
        'Rscript {params.script} '
        '--type gene_tpm '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_gene_vst:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_vst.tsv')
    shell:
        'Rscript {params.script} '
        '--type gene_vst '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '

rule tximport_gene_rlog:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_rlog.tsv')
    shell:
        'Rscript {params.script} '
        '--type gene_rlog '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_transcript_counts:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'transcript_counts.tsv')
    shell:
        'Rscript {params.script} '
        '--type tx '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_transcript_tpm:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'transcript_tpm.tsv')
    shell:
        'Rscript {params.script} '
        '--type tx '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '      

rule tximport_transcript_vst:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'transcript_vst.tsv')
    shell:
        'Rscript {params.script} '
        '--type tx_vst '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '

rule tximport_transcript_rlog:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'transcript_rlog.tsv')
    shell:
        'Rscript {params.script} '
        '--type tx_rlog '
        #'--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '

rule tximport_gene_info:
    input:
        rds = join(QUANT_INTERIM, '{quant}', 'tximport', SALMON_INDEX_TYPE + '_{quant}.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv'),
        gene_info = join(REF_DIR, 'anno', 'genes.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, '{quant}', 'tximport', 'gene_info.tsv')
    shell:
        'Rscript {params.script} '
        '--type gene_info '
        #'--txinfo {input.txinfo} '
        '--geneinfo {input.gene_info} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_transcript_info:
    input:
        rds = join(QUANT_INTERIM, 'salmon', 'tximport', SALMON_INDEX_TYPE + '_salmon.rds'),
        txinfo = join(REF_DIR, 'anno', 'transcripts.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, 'salmon', 'tximport', 'transcript_info.tsv')
    shell:
        'Rscript {params.script} '
        '--type tx_info '
        '--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
        
rule tximport_terminus_info:
    input:
        rds = join(QUANT_INTERIM, 'terminus', 'tximport', SALMON_INDEX_TYPE + '_terminus.rds'),
        txinfo = join(QUANT_INTERIM, 'terminus', 'terminus_info.tsv')
    params:
        script = src_gcf('scripts/tximport2csv.R')
    container:
        'docker://' + config['docker']['tximport']
    threads:
        8
    output:
        join(QUANT_INTERIM, 'terminus', 'tximport', 'transcript_info.tsv')
    shell:
        'Rscript {params.script} '
        '--type tx_info '
        '--txinfo {input.txinfo} '
        '--output {output} '
        '{input.rds} '
