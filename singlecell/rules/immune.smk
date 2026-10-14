# Parent-linked TCR uses external WT inputs and never includes expression rules.
import json
from pathlib import Path

TCR_ROOT = join(INTERIM_DIR, 'singlecell', 'immune')
TCR_ORDER = config['immune']['sublibrary_order']
TCR_PARENT_ROOT = config['immune']['wt_parent']['root']
TCR_SCRIPT = src_gcf('immune/scripts/tcr_workflow.py')
TCR_IMAGE = config['docker']['splitpipe_immune']
TCR_COMMIT = config.get('workflow_commit', 'UNRECORDED')
TCR_SETTINGS = json.dumps(config)
TCR_COMBINED_PARENT = join(TCR_PARENT_ROOT, config['immune']['wt_parent']['combined'])
TCR_COMBINED = join(TCR_ROOT, 'combined')


def tcr_parent(sublib):
    return join(TCR_PARENT_ROOT, config['immune']['wt_parent']['sublibraries'][sublib])


def tcr_parent_files(wc):
    # Declare metadata that Split-Pipe consumes, including per-sample cell calls.
    paths = []
    for parent in [tcr_parent(name) for name in TCR_ORDER] + [TCR_COMBINED_PARENT]:
        definition = join(parent, 'process', 'run_proc_def.json')
        paths.append(definition)
        if Path(definition).is_file():
            proc = json.loads(Path(definition).read_text())
            paths.extend(join(parent, sample['name'], 'DGE_filtered', 'cell_metadata.csv') for sample in proc['sample']['samples'])
        else:
            paths.append(join(parent, 'all-sample', 'DGE_filtered', 'cell_metadata.csv'))
    return paths


rule tcr_all:
    input:
        join(TCR_ROOT, 'linked', 'wt_tcr_cells.tsv.gz'),
        join(TCR_ROOT, 'linked', 'overlap.json'),
        join(TCR_COMBINED, 'complete.json'),
        join(TCR_COMBINED, 'all-sample_analysis_summary.html'),
        join(TCR_COMBINED, 'all-sample', 'TCR_unfiltered', 'tcr_annotation_airr.tsv.gz'),
        join(TCR_COMBINED, 'all-sample', 'TCR_filtered', 'barcode_report.tsv'),
        join(TCR_COMBINED, 'all-sample', 'TCR_filtered', 'clonotype_frequency.tsv'),
        join(TCR_COMBINED, 'all-sample', 'TCR_filtered', 'tcr_contigs.fa.gz')
    default_target: True


rule bfq_all:
    input:
        rules.tcr_all.input


rule tcr_preflight_all:
    input:
        expand(join(TCR_ROOT, 'preflight', '{sublib}', 'complete.json'), sublib=TCR_ORDER)


rule splitpipe_tcr_validate:
    input:
        tcr_parent_files
    output:
        manifest=join(TCR_ROOT, 'checks', 'wt_parents.json')
    params:
        script=TCR_SCRIPT,
        settings=TCR_SETTINGS
    threads: 1
    resources:
        mem_mb=2048
    log:
        join(TCR_ROOT, 'logs', 'validate_parents.log')
    container:
        'docker://' + TCR_IMAGE
    shell:
        'python {params.script:q} validate-parents --settings {params.settings:q} --output {output.manifest:q} > {log:q} 2>&1'


rule splitpipe_tcr_preflight:
    input:
        unpack(get_raw_fastq),
        parents=rules.splitpipe_tcr_validate.output.manifest
    output:
        complete=join(TCR_ROOT, 'preflight', '{sublib}', 'complete.json')
    params:
        script=TCR_SCRIPT,
        output=lambda wc: str(Path(join(TCR_ROOT, 'preflight', wc.sublib)).resolve()),
        parent=lambda wc: tcr_parent(wc.sublib),
        reference=config['immune']['reference'],
        provenance=lambda wc: join(TCR_ROOT, 'provenance', 'preflight_' + wc.sublib + '.json'),
        image=TCR_IMAGE,
        commit=TCR_COMMIT
    threads: 12
    resources:
        mem_mb=96000,
        io_jobs=1
    log:
        join(TCR_ROOT, 'logs', 'preflight_{sublib}.log')
    benchmark:
        join(TCR_ROOT, 'benchmarks', 'preflight_{sublib}.tsv')
    container:
        'docker://' + TCR_IMAGE
    shell:
        'python {params.script:q} run --dryrun --name {wildcards.sublib:q} --fastqs {input.R1:q} {input.R2:q} '
        '--parent {params.parent:q} --reference {params.reference:q} --output-dir {params.output:q} '
        '--threads {threads} --image {params.image:q} --workflow-commit {params.commit:q} '
        '--provenance {params.provenance:q} --success {output.complete:q} > {log:q} 2>&1'


rule splitpipe_tcr:
    input:
        unpack(get_raw_fastq),
        parents=rules.splitpipe_tcr_validate.output.manifest,
        preflight=join(TCR_ROOT, 'preflight', '{sublib}', 'complete.json')
    output:
        complete=join(TCR_ROOT, 'sublibraries', '{sublib}', 'complete.json'),
        runproc=join(TCR_ROOT, 'sublibraries', '{sublib}', 'process', 'run_proc_def.json'),
        airr=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_filtered', 'tcr_annotation_airr.tsv.gz'),
        raw_airr=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_unfiltered', 'tcr_annotation_airr.tsv.gz'),
        report=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample_analysis_summary.html'),
        barcode=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_filtered', 'barcode_report.tsv'),
        raw_barcode=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_unfiltered', 'barcode_report.tsv'),
        contigs=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_filtered', 'tcr_contigs.fa.gz'),
        raw_contigs=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_unfiltered', 'tcr_contigs.fa.gz'),
        clones=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_filtered', 'clonotype_frequency.tsv'),
        raw_clones=join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_unfiltered', 'clonotype_frequency.tsv')
    params:
        script=TCR_SCRIPT,
        output=lambda wc: str(Path(join(TCR_ROOT, 'sublibraries', wc.sublib)).resolve()),
        parent=lambda wc: tcr_parent(wc.sublib),
        reference=config['immune']['reference'],
        provenance=lambda wc: join(TCR_ROOT, 'provenance', 'run_' + wc.sublib + '.json'),
        image=TCR_IMAGE,
        commit=TCR_COMMIT
    threads: 12
    resources:
        mem_mb=96000,
        io_jobs=1
    log:
        join(TCR_ROOT, 'logs', 'run_{sublib}.log')
    benchmark:
        join(TCR_ROOT, 'benchmarks', 'run_{sublib}.tsv')
    container:
        'docker://' + TCR_IMAGE
    shell:
        'python {params.script:q} run --name {wildcards.sublib:q} --fastqs {input.R1:q} {input.R2:q} '
        '--parent {params.parent:q} --reference {params.reference:q} --output-dir {params.output:q} '
        '--threads {threads} --image {params.image:q} --workflow-commit {params.commit:q} '
        '--provenance {params.provenance:q} --success {output.complete:q} > {log:q} 2>&1'


rule splitpipe_tcr_combine:
    input:
        complete=expand(join(TCR_ROOT, 'sublibraries', '{sublib}', 'complete.json'), sublib=TCR_ORDER),
        definitions=expand(join(TCR_ROOT, 'sublibraries', '{sublib}', 'process', 'run_proc_def.json'), sublib=TCR_ORDER),
        filtered=expand(join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_filtered', 'tcr_annotation_airr.tsv.gz'), sublib=TCR_ORDER),
        unfiltered=expand(join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_unfiltered', 'tcr_annotation_airr.tsv.gz'), sublib=TCR_ORDER),
        supporting=expand(join(TCR_ROOT, 'sublibraries', '{sublib}', 'all-sample', 'TCR_{fraction}', '{filename}'),
                          sublib=TCR_ORDER, fraction=['filtered', 'unfiltered'],
                          filename=['barcode_report.tsv', 'tcr_contigs.fa.gz', 'clonotype_frequency.tsv']),
        parents=rules.splitpipe_tcr_validate.output.manifest
    output:
        complete=join(TCR_COMBINED, 'complete.json'),
        airr=join(TCR_COMBINED, 'all-sample', 'TCR_filtered', 'tcr_annotation_airr.tsv.gz'),
        report=join(TCR_COMBINED, 'all-sample_analysis_summary.html'),
        raw_airr=join(TCR_COMBINED, 'all-sample', 'TCR_unfiltered', 'tcr_annotation_airr.tsv.gz'),
        barcode=join(TCR_COMBINED, 'all-sample', 'TCR_filtered', 'barcode_report.tsv'),
        contigs=join(TCR_COMBINED, 'all-sample', 'TCR_filtered', 'tcr_contigs.fa.gz'),
        clones=join(TCR_COMBINED, 'all-sample', 'TCR_filtered', 'clonotype_frequency.tsv')
    params:
        script=TCR_SCRIPT,
        output=str(Path(TCR_COMBINED).resolve()),
        sublibraries=[str(Path(join(TCR_ROOT, 'sublibraries', name)).resolve()) for name in TCR_ORDER],
        parent=TCR_COMBINED_PARENT,
        provenance=join(TCR_ROOT, 'provenance', 'combine.json'),
        image=TCR_IMAGE,
        commit=TCR_COMMIT
    threads: 12
    resources:
        mem_mb=96000,
        io_jobs=1
    log:
        join(TCR_ROOT, 'logs', 'combine.log')
    benchmark:
        join(TCR_ROOT, 'benchmarks', 'combine.tsv')
    container:
        'docker://' + TCR_IMAGE
    shell:
        'python {params.script:q} run --sublibraries {params.sublibraries:q} --parent {params.parent:q} '
        '--output-dir {params.output:q} --threads {threads} --image {params.image:q} --workflow-commit {params.commit:q} '
        '--provenance {params.provenance:q} --success {output.complete:q} > {log:q} 2>&1'


rule splitpipe_tcr_link_wt:
    input:
        airr=rules.splitpipe_tcr_combine.output.airr,
        complete=rules.splitpipe_tcr_combine.output.complete,
        cells=join(TCR_COMBINED_PARENT, 'all-sample', 'DGE_filtered', 'cell_metadata.csv')
    output:
        table=join(TCR_ROOT, 'linked', 'wt_tcr_cells.tsv.gz'),
        summary=join(TCR_ROOT, 'linked', 'overlap.json')
    params:
        script=TCR_SCRIPT
    threads: 1
    resources:
        mem_mb=16000
    log:
        join(TCR_ROOT, 'logs', 'link_wt.log')
    container:
        'docker://' + TCR_IMAGE
    shell:
        'python {params.script:q} link-wt --airr {input.airr:q} --metadata {input.cells:q} '
        '--output {output.table:q} --summary {output.summary:q} > {log:q} 2>&1'
