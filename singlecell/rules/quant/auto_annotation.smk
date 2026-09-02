#-*- mode:snakemake -*-
"""
Automatic celltype  annotation of single cell rna-seq data
"""
include:
    join(GCFDB_DIR, 'allen_institute.smk')
include:
    join(GCFDB_DIR, 'celltypist.smk')


MM_ORG = config.get('celltype_annotation', {}).get('orthologs')
MM_ORG = MM_ORG or config['organism']

PRE_ANNO = True
_AGGR_ID = config['quant']['aggregate']['groupby']

rule orthogene_premap:
    input:
        unpack(get_raw_mtx)
    output:
        gene_map = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation',  'orthogene', 'orthologs.tsv')
    params:
        script = src_gcf('scripts/run_orthogene.R'),
        src_org = config['organism'],
        dst_org = MM_ORG,
        method = 'gprofiler',
        non121_strategy = 'drop_both_species'
    container:
        'docker://gcfntnu/orthogene:1.12.0'
    threads:
        24
    shell:
        'Rscript {params.script} '
        '--input {input.mtx} '
        '--output {output.gene_map} '
        '--src {params.src_org} '
        '--dst {params.dst_org} '
        '--no-cache '


rule orthogene_premap_aggr:
    input:
        expand(join(QUANT_INTERIM, '{{quantifier}}', '{sample}', 'annotation',  'orthogene', 'orthologs.tsv'), sample=AGGR_IDS[_AGGR_ID])
    output:
        tsv = join(QUANT_INTERIM, 'aggregate', '{quantifier}', '{aggr_id}_orthologs.tsv')
    params:
        script = src_gcf('scripts/aggr_orthogene.py')
    container:
        'docker://' + config['docker']['default']
    threads:
        24
    shell:
        'python {params.script} '
        '--output {output.tsv} '
        ' {input} '
        
rule mapmycells_clean_premap_input:
    input:
        unpack(get_filtered_mtx),
        gene_map = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation', 'orthogene', 'orthologs.tsv')
    output:
        tiny_h5ad = temp(join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation', 'mapmycells', '_tiny.h5ad'))
    params:
        script = src_gcf('scripts/mapmycells_input.py'),
        src_organism = config['organism'],
        dst_organism = 'mus_musculus' if config['organism'] in ['mus_musculus', 'rattus_norvegicus'] else 'homo_sapiens'
    container:
        'docker://' + config['docker']['scanpy']
    threads:
        80
    shell: 
        'python {params.script} '
        '--input {input.mtx} '
        '--output {output.tiny_h5ad} '
        '--gene-map {input.gene_map} '
        '--src-organism {params.src_organism} '
        '--dst-organism {params.dst_organism}'
    

rule mapmycells_premap_from_specified_markers:
    input:
        tiny_h5ad = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation',  'mapmycells', '_tiny.h5ad'),
        pre_stats_h5 = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', MM_ORG, 'precomputed_stats.h5'),
        markers_json = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', MM_ORG, 'markers.json')
    output:
        anno_csv = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation',  'mapmycells', 'annotation.csv'),
        anno_json = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation',  'mapmycells', 'annotation.json')
    params:
        args = '--type_assignment.chunk_size 3000 --type_assignment.bootstrap_factor 0.5 --type_assignment.bootstrap_iteration 100 --type_assignment.normalization raw --type_assignment.rng_seed 661123 '
    container:
        'docker://gcfntnu/mapmycells:1.5.1'
    threads:
        48
    shell:
        'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; '
        'python -m cell_type_mapper.cli.from_specified_markers '
        '--precomputed_stats.path {input.pre_stats_h5} '
        '--query_markers.serialized_lookup {input.markers_json} '
        '--type_assignment.n_processors {threads} '
        '--query_path {input.tiny_h5ad} '
        '--extended_result_path {output.anno_json} '
        '--csv_result_path {output.anno_csv} '
        '--tmp_dir /dev/shm/mapmycells_{wildcards.sample} '
        ' {params.args} '

    
rule mapmycells_output_processing:
    input:
        anno_csv = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation',  'mapmycells', 'annotation.csv'),
        colors = join(EXT_DIR, 'allen-brain-cell-atlas', 'derived', 'abc_colors.json'),
        meta = 'data/ext/cl.df_CCN202307220.xlsx'
    output:
        extended_anno_csv = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'annotation',  'mapmycells', 'annotation_extended.csv')
    params:
        script = src_gcf("scripts/mapmycells_colormap.py")
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '--annotation {input.anno_csv} '
        '--colors {input.colors} '
        '--metadata {input.meta} '
        '--preset minimal '
        '--out {output.extended_anno_csv} '
        '--verbose '

def aggr_input(wildcards):
    samples_by_aggr_id = AGGR_IDS.get(wildcards.aggr_id)
    if config.get('celltype_annotation', {}).get('mapmycells', {}).get('extended', False):
        input_files = expand(rules.mapmycells_output_processing.output.extended_anno_csv,
                             quantifier=wildcards.method,
                             sample=samples_by_aggr_id)
    else:
        input_files = expand(rules.mapmycells_premap_from_specified_markers.output.anno_csv,
                             quantifier=wildcards.method,
                             sample=samples_by_aggr_id)
    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        aggr_csv = join(QUANT_INTERIM, 'aggregate', 'description', f'{wildcards.aggr_id}_aggr.csv')
        return {'input_files': input_files, 'aggr_csv': aggr_csv}
    else:
        return {'input_files': input_files}


rule mapmycells_premap_aggr:
    input:
        unpack(aggr_input),
    output:
        join(QUANT_INTERIM, 'aggregate', '{method}' , '{aggr_id}_premap_annotation.tsv')
    params:
        script = src_gcf("scripts/aggr_barcode_info.py"),
        args = barcode_aggr_args
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '{input.input_files} '
        '{params.args} '
        '--output {output} '
        
        
rule mapmycells_gene_map:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_raw.h5ad')
    output:
        gene_map = join(QUANT_INTERIM, 'aggregate', '{method}', '{aggr_id}_gene_map.tsv')
    params:
        script = src_gcf('scripts/run_orthogene.R'),
        src_org = config['organism'],
        dst_org = MM_ORG,
        method = 'gprofiler',
        non121_strategy = 'drop_both_species'
    container:
        'docker://gcfntnu/orthogene:1.12.0'
    shell:
        'Rscript {params.script} '
        '{input.aggr_raw_h5ad} '
        '{output.gene_map} '
        '{params.src_org} '
        '{params.dst_org} '

        
rule mapmycells_clean_input:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad'),
        gene_map = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.gene_map.tsv')
    output:
        aggr_clean_h5ad = temp(join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_tiny.h5ad'))
    params:
        script = src_gcf('scripts/mapmycells_input.py'),
        src_organism = config['organism'],
        dst_organism = 'mus_musculus' if config['organism'] in ['mus_musculus', 'rattus_norvegicus'] else 'homo_sapiens'
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} '
        '--input {input.aggr_raw_h5ad} '
        '--output {output.aggr_clean_h5ad} '
        '--gene-map {input.gene_map} '
        '--src-organism {params.src_organism} '
        '--dst-organism {params.dst_organism} '
    

rule mapmycells_from_specified_markers:
    input:
        aggr_filtered_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_tiny.h5ad'),
        pre_stats_h5 = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', MM_ORG, 'precomputed_stats.h5'),
        markers_json = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', MM_ORG, 'markers.json')
    output:
        anno_csv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_annotate', '{aggr_id}_mapmycells_annotation.csv'),
        anno_json = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_annotate', '{aggr_id}_mapmycells_annotation.json')
    params:
        args = '--type_assignment.bootstrap_factor 0.5 --type_assignment.bootstrap_iteration 100 --type_assignment.normalization raw --type_assignment.rng_seed 661123 '
    container:
        'docker://gcfntnu/mapmycells:1.5.1'
    threads:
        64
    shell:
        'python -m cell_type_mapper.cli.from_specified_markers '
        '--precomputed_stats.path {input.pre_stats_h5} '
        '--query_markers.serialized_lookup {input.markers_json} '
        '--type_assignment.n_processors {threads} '
        '--query_path {input.aggr_raw_h5ad} '
        '--extended_result_path {output.anno_json} '
        '--csv_result_path {output.anno_csv} '
        '{params.args} '


rule celltypist_default_models:
    params:
        script = src_gcf('scripts/download_celltypist_models.py'),
        celltypist_folder = join(EXT_DIR, 'celltypist')
    output:
        model = join(EXT_DIR, 'celltypist', 'data', 'models', 'Healthy_COVID19_PBMC.pkl')
    container:
        'docker://gcfntnu/rapids-scanpy:latest'
    shell:
        'export CELLTYPIST_FOLDER="{params.celltypist_folder}" '
        '&& '
        'python -c "from celltypist import models; models.download_models(force_update=True)"'


rule run_celltypist:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad'),
        model = join(EXT_DIR, 'celltypist', 'data', 'models', config.get('celltype_annotation', {}).get('celltypist', {}).get('model', '')),
        qc_mask = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_qc', '{aggr_id}_sctk_autoqc_mask.tsv'),
        gene_map = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.gene_map.tsv')
    output:
        anno_csv = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_annotate', '{aggr_id}_celltypist_annotation.csv')
    params:
        script = src_gcf('scripts/run_celltypist.py'),
        src_organism = config['organism'],
        dst_organism = 'mus_musculus' if config['organism'] in ['mus_musculus', 'rattus_norvegicus'] else 'homo_sapiens',
        args = '--use-GPU --no-qc-filter --plot '
    container:
        'docker://gcfntnu/rapids-scanpy:latest'
    shell:
        'python {params.script} '
        '--input {input.aggr_raw_h5ad} '
        '--model {input.model} '
        '--output {output.anno_csv} '
        '--gene-map {input.gene_map} '
        '--qc-mask {input.qc_mask} '
        '--src-organism {params.src_organism} '
        '--dst-organism {params.dst_organism} '
        '{params.args} '
    

rule auto_annotate_scanpy:
    input:
        aggr_raw_h5ad = join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_filtered.h5ad'),
        anno = expand(join(QUANT_INTERIM, 'aggregate', '{{method}}', 'auto_annotate', '{{aggr_id}}_{anno_method}_annotation.csv'), anno_method=config['celltype_annotation']['method'].split(','))
    output:
        aggr_anno = join(QUANT_INTERIM, 'aggregate', '{method}', 'auto_annotate', '{aggr_id}_celltype_annotation.tsv')
    params:
        script = src_gcf('scripts/combine_celltype_annotations.py')
    container:
        'docker://' + config['docker']['scanpy']
    log:
        'logs/{method}/{aggr_id}_auto_annotate_scanpy.log'
    shell:
        'python {params.script} '
        '--input {input.aggr_raw_h5ad} '
        '--output {output.aggr_anno} '
        '--log {log} '
        '--annotation {input.anno} '

