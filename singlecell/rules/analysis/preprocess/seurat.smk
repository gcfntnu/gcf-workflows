

rule seurat_filter_cells:
    container:
        'docker://quay.io/biocontainers/r-seurat-scripts:0.0.3--r341_0'
        
rule seurat_filter_genes:
    container:
        'docker://quay.io/biocontainers/r-seurat-scripts:0.0.3--r341_0'

rule seurat_normalize:
    container:
        'docker://quay.io/biocontainers/r-seurat-scripts:0.0.3--r341_0'

rule seurat_variable_genes:
    container:
        'docker://quay.io/biocontainers/r-seurat-scripts:0.0.3--r341_0'
