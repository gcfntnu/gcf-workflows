#-*- mode:snakemake -*-

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

include:
    join(GCFDB_DIR, 'unitas.db')


rule unitas_mirtrace_fasta:
    input:
        join(FILTER_INTERIM, 'mirtrace', 'qc_passed_reads.all.collapsed', '{sample}_R1.fasta')
    output:
        join(FILTER_INTERIM, 'mirtrace', 'unitas_fasta', '{sample}_R1.fasta')
    shell:
        """
        sed -r '/>/ s/(.*)x([[:digit:]]+)(.*)/>\\2/' {input} > {output}
        """

def get_unitas_seq(wildcards):
    if config['filter']['trim']['quantifier'] == 'mirtrace':
        fastq = join(FILTER_INTERIM, 'mirtrace', 'unitas_fasta', '{sample}_R1.fasta')
    else:
        fastq = get_filtered_fastq(wildcards)['R1']
    return fastq

def get_spikein_fasta():
    ref = config['filter']['spikein'].get('ref', '').lower()
    if ref == 'ercc':
        return join(EXT_DIR, 'ERCC', 'fasta', 'ERCC92.fa')
    elif ref == 'smallrna_calibrators':
        return join(EXT_DIR, 'spikein', 'fasta', 'smallrna_calibrators.fa')
    else:
        return join(EXT_DIR, 'spikein', 'fasta', 'merged_spikein.fa')
    
rule spikein_unitas:
    input:
        get_spikein_fasta()
    params:
        script = src_gcf('scripts/spikein_unitas.py')
    output:
        join(QUANT_INTERIM, 'unitas', '_spikein_unitas_formatted.fa')
    shell:
        'python {params.script} {input} > {output}'
        
UNITAS_REFSEQ = []
if config['filter']['spikein']['quantifier'] == 'unitas':
    UNITAS_REFSEQ.append(rules.spikein_unitas.output)
if config['filter']['contaminants']['quantifier'] == 'unitas':
    UNITAS_REFSEQ.append(rules.contaminants_unitas.output)
    
rule unitas_extra_reference:
    input:
        UNITAS_REFSEQ
    output:
        join(QUANT_INTERIM, 'unitas', 'refseq.fa')
    run:
        if len(input) > 0:
            shell('cat {input} > {output}')
        else:
            shell('touch {output}')
    
UNITAS_ARGS = ''
if config['filter']['ribosomal']['quantifier'] == 'unitas':
    UNITAS_ARGS += '--riborase '

if len(UNITAS_REFSEQ) > 0:
    UNITAS_ARGS += '-refseq {} '.format(abspath(join(QUANT_INTERIM, 'unitas', 'refseq.fa')))


rule unitas_rundir:
    input:
        seqmap = rules.unitas_seqmap.output,
        db_log = rules.unitas_db.log,
        db = rules.unitas_db.output.db,
    output:
        seqmap = 'seqmap.exe'
    params:
        outdir = rules.unitas_db.output.db,
        seqmap = temp('seqmap.exe'),
        db = directory('UNITAS_refdump_{}'.format(ORG)),
    shell:
        """
        cp {input.seqmap} .
        ln -srf {params.outdir} .
        """
        
rule unitas:
    input:
        seqmap = 'seqmap.exe',
        fastq = get_unitas_seq,
        refseq = rules.unitas_extra_reference.output
    params:
        outdir = abspath(join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}')),
        args = UNITAS_ARGS,
        org  = config['organism']
    output:
        html = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'results.html'),
        trftable = join(QUANT_INTERIM, 'unitas',  'unitas', '{sample}', 'unitas.tRF-table.simplified.txt'),
        annotation = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'unitas.annotation_summary.txt'),
        full_annotation = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'unitas.full_annotation_matrix.txt'),
        target_hits = join(QUANT_INTERIM, 'unitas',  'unitas', '{sample}', 'unitas.hits_per_target.txt'),
        modifications = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'unitas.miR-modifications_{}.txt'.format(UNITAS_ORG)),
        trftable_ext = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'unitas.tRF-table.txt'),
        mirtable = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'miR-table_{}.simplified.txt'.format(UNITAS_ORG)),
        pirna = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'unitas.piRcandidates-tails_{}.txt'.format(UNITAS_ORG)),
        mirtable_ext= join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}', 'unitas.miR-table_{}.txt'.format(UNITAS_ORG))
    container:
        'docker://' + config['docker']['unitas']
    threads:
        2
    shell:
        """
        unitas.pl -input {input.fastq} -s {params.org} -threads {threads} {params.args} 
        rm -rf {params.outdir}; mkdir -p {params.outdir}
        mv -f UNITAS_*{wildcards.sample}_R1.fast*_#1/* {params.outdir}/
        rm -rf UNITAS_*{wildcards.sample}_R1.*fas*_#1
        """

rule unitas_annotation2yaml:
    input:
        rules.unitas.output.annotation
    params:
        script = src_gcf('scripts/unitas_annotation2json.py')
    output:
        join(QUANT_INTERIM, 'unitas', '{sample}', 'annotation_summary.json')
    shell:
        'python {params.script} {input} > {output}'

rule unitas_mirtable:
    input:
        expand(rules.unitas.output.mirtable, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_mirtable.py')
    output:
        counts = join(QUANT_INTERIM, 'unitas', 'mir_counts.tsv')
    shell:
        'python {params.script} {input} > {output}'
        
rule unitas_isomirtable:
    input:
        expand(rules.unitas.output.mirtable_ext, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_isomirtable.py'),
        min_count = 0
    output:
        counts = join(QUANT_INTERIM, 'unitas', 'isomir_counts.tsv'),
        anno = join(QUANT_INTERIM, 'unitas', 'isomir_counts_anno.tsv')
    shell:
        'python {params.script} {input} --min-count {params.min_count} --output {output.counts}'

rule unitas_trftable:
    input:
        expand(rules.unitas.output.trftable, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_trftable.py')
    output:
        join(QUANT_INTERIM, 'unitas', 'trf_counts.tsv')
    shell:
        'python {params.script} {input} > {output}'
        
rule unitas_isotrftable:
    input:
        expand(rules.unitas.output.trftable_ext, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_isotrftable.py')
    output:
        join(QUANT_INTERIM, 'unitas', 'isotrf_counts.csv')
    shell:
        'python {params.script} {input} > {output}'

rule unitas_allhits:
    input:
        expand(rules.unitas.output.target_hits, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_allhits.py') 
    output:
        join(QUANT_INTERIM, 'unitas', 'allhits.tsv')
    shell:
        'python {params.script} {input} --output {output} '
        
rule unitas_modifications:
    input:
        expand(rules.unitas.output.modifications, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_modifications.py'),
        outdir = join(QUANT_INTERIM, 'unitas')
    output:
        join(QUANT_INTERIM, 'unitas', 'modification_per_position.tsv'),
        join(QUANT_INTERIM, 'unitas', 'internal_modification.tsv'),
        join(QUANT_INTERIM, 'unitas', 'internal_modification_position.tsv'),
        join(QUANT_INTERIM, 'unitas', 'NT_3p_extension.tsv')
    shell:
        'python {params.script} -o {params.outdir} {input}'
        
rule unitas_annotations:
    input:
        expand(rules.unitas.output.annotation, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_annotations.py')
    output:
        join(QUANT_INTERIM, 'unitas', 'annotations.tsv')
    shell:
        'python {params.script} {input} > {output}'

rule unitas_allfeatures:
    input:
        expand(rules.unitas.output.full_annotation, sample=SAMPLES)
    params:
        script = src_gcf('scripts/unitas_feature_counts.py') 
    output:
        join(QUANT_INTERIM, 'unitas', 'allfeatures.tsv')
    shell:
        'python {params.script} --output {output} {input} '

rule unitas_isomir_tximport:
    input:
        counts = rules.unitas_isomirtable.output.counts,
        anno = rules.unitas_isomirtable.output.anno
    output:
       join(QUANT_INTERIM, 'unitas', 'tximport.rds')
    params:
        script = src_gcf('scripts/isomir2tximport.R')
    shell:
        'Rscript {params.script} '
        '--input {input.counts} '
        '--feature-info {input.anno} '
        '--type unitas '
        '--output {output} -v '

rule unitas_mirna_anndata:
    input:
        mir_table = rules.unitas_mirtable.output,
        sample_info = join(INTERIM_DIR, 'sample_info.tsv'),
        feature_info = rules. mirbase_genes.output
    output:
        join(QUANT_INTERIM, 'unitas', 'adata.h5ad')
    params:
        script = src_gcf('scripts/create_anndata.py')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--input {input.mir_table} '
        '--sample-info {input.sample_info} '
        '--feature-info {input.feature_info} '
        '--add-vsn '
        '--output {output} '

rule unitas_biotypes:
    input:
        
rule unitas_all:
    input:
        rules.unitas_mirtable.output,
        rules.unitas_trftable.output,
        rules.unitas_isomirtable.output,
        rules.unitas_isotrftable.output,
        rules.unitas_modifications.output,
        rules.unitas_annotations.output

