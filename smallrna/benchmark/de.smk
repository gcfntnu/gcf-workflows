#-*- mode:snakemake -*-
"""Differential expression on a ratiometric pool's evaluation (#414): de.R in the `sartools` image of docker.config, for
the whole pool and for each level of `within` (#444) that `benchmark.py evaluate` listed in benchmark/de/levels.tsv.

Run in the evaluated project after `benchmark.py evaluate` has written benchmark/de/, with the same container settings
as the project's `quant_all` (the repository holding de.R must be visible in the container):

    snakemake -s <repo>/smallrna/benchmark/de.smk -d <project> --use-apptainer --cores 1

then score the results with `benchmark.py evaluate --project <project> --de`.
"""
from pathlib import Path

import yaml
from snakemake.exceptions import WorkflowError

HERE = Path(str(workflow.current_basedir))
SARTOOLS = yaml.safe_load((HERE.parents[1] / 'docker.config').read_text())['docker']['sartools']
DE = 'benchmark/de'
DE_INPUT = DE + '/input'
LEVELS_TABLE = Path(DE, 'levels.tsv')
if not LEVELS_TABLE.is_file():
    raise WorkflowError('{} not found in {}: run `benchmark.py evaluate` in the project first; it writes the inputs of '
                        'differential expression'.format(LEVELS_TABLE, Path.cwd()))
LEVELS = [line.split('\t') for line in LEVELS_TABLE.read_text().splitlines()[1:] if line]  # [column, level]
COMMAND = ('Rscript {params.script:q} --counts {input.counts:q} --library-sizes {input.library_sizes:q} --design {input.design:q} '
           '--results {output.results:q} --record {output.record:q} > {log:q} 2>&1')


rule all:
    input:
        DE + '/results.tsv',
        ['{}/{}/{}/results.tsv'.format(DE, column, level) for column, level in LEVELS]


rule benchmark_de:
    input:
        counts = DE_INPUT + '/counts.tsv',
        library_sizes = DE_INPUT + '/library_sizes.tsv',
        design = DE_INPUT + '/design.tsv',
        script = str(HERE / 'de.R')  # an edit to de.R reruns the fits, as evaluate --de checks its MD5
    output:
        results = DE + '/results.tsv',
        record = DE + '/record.tsv'
    params:
        script = str(HERE / 'de.R')
    threads:
        1
    resources:
        mem_mb = 2000
    log:
        DE + '/de.log'
    container:
        'docker://' + SARTOOLS
    shell:
        COMMAND


rule benchmark_de_level:
    input:
        counts = DE + '/{column}/{level}/input/counts.tsv',
        library_sizes = DE + '/{column}/{level}/input/library_sizes.tsv',
        design = DE + '/{column}/{level}/input/design.tsv',
        script = str(HERE / 'de.R')
    output:
        results = DE + '/{column}/{level}/results.tsv',
        record = DE + '/{column}/{level}/record.tsv'
    params:
        script = str(HERE / 'de.R')
    threads:
        1
    resources:
        mem_mb = 2000
    log:
        DE + '/{column}/{level}/de.log'
    wildcard_constraints:
        column = '[^/]+',
        level = '[^/]+'
    container:
        'docker://' + SARTOOLS
    shell:
        COMMAND
