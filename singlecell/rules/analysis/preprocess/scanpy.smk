

rule scanpy_concat:
    input:
        expand(join(QUANT_INTERIM, 'star', '{sample}', 'Solo.out', 'matrix.mtx'), sample=SAMPLES)
    output:
        join(ANALYSIS_INTERIM, 'preprocess', 'scanpy', 'data.5ahd')
    params:
        script = srcdir('scripts/scanpy_read.py')
    shell:
        'python {params.script} '
        '-i {input} '
        '-o {output} '
        
rule scanpy_filter_cells:
    input:
        rules.scanpy_concat.output
    output:
        join(ANALYSIS_INTERIM, 'preprocess', 'scanpy', 'filtered1.5ahd')
    singularity:
        'quay.io/biocontainers/scanpy-scripts:0.0.4--py37_1'
    shell:
        'scanpy-filter-cells.py '
        '-i {input} '
        '-f auto-detect '
        '-F anndata '
        '-l 300 '
           
rule scanpy_filter_genes:
    input:
        rules.scanpy_filter_cells.output
    output:
        join(ANALYSIS_INTERIM, 'preprocess', 'scanpy', 'filtered2.5ahd')
    singularity:
        'quay.io/biocontainers/scanpy-scripts:0.0.4--py37_1'
    shell:
        'scanpy-filter-genes.py '
        '-i {input} '
        
rule scanpy_normalize:
    singularity:
        'quay.io/biocontainers/scanpy-scripts:0.0.4--py37_1'

rule scanpy_variable_genes:
    singularity:
        'quay.io/biocontainers/scanpy-scripts:0.0.4--py37_1'

