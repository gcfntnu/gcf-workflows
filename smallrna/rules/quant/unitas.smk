#-*- mode:snakemake -*-

include:
    join(GCFDB_DIR, 'unitas.db')


def get_unitas_seq(wildcards):
    return get_filtered_fastq(wildcards)['R1']

SPIKEIN_SETS = ('ERCC92', 'smallrna_calibrators')

def get_spikein_fasta():
    """FASTA files of the comma-separated spike-in sets in filter.spikein.ref."""
    ref = config['filter']['spikein'].get('ref')
    names = [name.strip() for name in str(ref or '').split(',') if name.strip()]
    unknown = [name for name in names if name not in SPIKEIN_SETS]
    if unknown:
        raise ValueError('filter.spikein.ref: unsupported value {!r}; use a comma-separated list of {}'.format(ref, ', '.join(SPIKEIN_SETS)))
    if not names and config['filter']['spikein']['quantifier'] == 'unitas':
        raise ValueError('filter.spikein.quantifier is unitas but filter.spikein.ref lists no spike-in sets')
    return [join(EXT_DIR, 'spikein', 'fasta', name + '.fa') for name in names]
    
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


rule unitas:
    input:
        fastq = get_unitas_seq,
        refseq = rules.unitas_extra_reference.output,
        seqmap = rules.unitas_seqmap.output,
        db = rules.unitas_db.output.db,
        db_log = rules.unitas_db.log
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
    params:
        script = src_gcf('scripts/run_unitas.sh'),
        workdir = join(QUANT_INTERIM, 'unitas', 'work', '{sample}'),
        outdir = join(QUANT_INTERIM, 'unitas', 'unitas', '{sample}'),
        org = config['organism'],
        args = UNITAS_ARGS
    threads:
        2
    container:
        'docker://' + config['docker']['unitas']
    shell:
        'bash {params.script} --input {input.fastq} --species {params.org} --threads {threads} '
        '--seqmap {input.seqmap} --refdump {input.db} '
        '--workdir {params.workdir} --outdir {params.outdir} -- {params.args}'

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

rule unitas_mirna_vsn:
    input:
        rules.unitas_mirtable.output
    output:
        join(QUANT_INTERIM, 'unitas', 'mirna_vsd.tsv')
    params:
        script = src_gcf('scripts/mirna_vsn.R')
    container:
        'docker://' + config['docker']['tximport']
    shell:
        'Rscript {params.script:q} --input {input:q} --output {output:q}'

rule unitas_mirna_anndata:
    input:
        mir_table = rules.unitas_mirtable.output,
        sample_info = join(INTERIM_DIR, 'sample_info.tsv'),
        feature_info = rules.mirbase_genes.output,
        vsn = rules.unitas_mirna_vsn.output
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
        '--add-vsn --vsn-input {input.vsn:q} '
        '--output {output} '

rule unitas_all:
    input:
        rules.unitas_mirtable.output,
        rules.unitas_trftable.output,
        rules.unitas_isomirtable.output,
        rules.unitas_isotrftable.output,
        rules.unitas_modifications.output,
        rules.unitas_annotations.output

