rule seurat_cellcycle:
    input:
        seurat = 'seurat.rds',
        markers = 'data/ext/regev_lab_cell_cycle_genes.txt'
    output:
        'seurat_cellcycle.rds'
    params:
        script = workflow.source_path('scripts/seurat_cellcycle.R'),
        args = '--write-scores -v'
    container:
        ''
    shell:
        'Rscript {params.script} '
        '-i {input.seurat} '
        '-m {input.markers} '
        '-o {output} '
        '{params.args} '
        
