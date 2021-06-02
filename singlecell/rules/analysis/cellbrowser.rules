
rule cellbrowser_cellranger:
    input:
        rules.cellranger_aggr.output.filt_mtx
    output:
        join(QUANT_INTERIM, 'cellbrowser', 'cellranger')
    params:
        name = config['project_id']
    singularity:
        'docker://' + config['docker']['cellbrowser']
    shell:
        'cbImportCellranger -i {input} -o {output} --name {params.name}'

rule cellbrowser_scanpy:
    input:
        rules.cellranger_aggr.output.filt_mtx
    output:
        join(QUANT_INTERIM, 'cellbrowser', 'scanpy')
    params:
        name = config['project_id']
    singularity:
        'docker://' + config['docker']['cellbrowser']
    shell:
        'cbImportCellranger -i {input} -o {output} --name {params.name}'

