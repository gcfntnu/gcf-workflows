#-*- mode:snakemake -*-
"""Differential expression on a benchmark evaluation (#414, #447): de.R in the `sartools` image of docker.config, for
the whole dataset's DE set and for each other set that `benchmark.py evaluate` listed in benchmark/de/levels.tsv (the
levels of a pool's `within`, #444; the synthetic dataset's set with the nuisance covariate; a biological dataset's null
comparisons, #449). A biological dataset that declares null comparisons only has no whole set: evaluate writes no
benchmark/de/input/.

Run in the evaluated project after `benchmark.py evaluate` has written benchmark/de/, with the same container settings
as the project's `quant_all` (the repository holding de.R must be visible in the container):

    snakemake -s <repo>/smallrna/benchmark/de.smk -d <project> --use-apptainer --cores 1

then score the results with `benchmark.py evaluate --project <project> --de`.

`--config de=<folder>` runs the sets of another DE folder, relative to the working directory: `benchmark.py depth`
writes those of each depth to de/ in its output (#414), run with `-d <depth output> --config de=de`, then scored with
`benchmark.py depth --de`.
"""
from pathlib import Path

import yaml
from snakemake.exceptions import WorkflowError

HERE = Path(str(workflow.current_basedir))
SARTOOLS = yaml.safe_load((HERE.parents[1] / 'docker.config').read_text())['docker']['sartools']
DE = str(config.get('de', 'benchmark/de'))
DE_INPUT = DE + '/input'
LEVELS_TABLE = Path(DE, 'levels.tsv')
if not LEVELS_TABLE.is_file():
    raise WorkflowError('{} not found in {}: run `benchmark.py evaluate` in the project first (`depth` for --config '
                        'de=de in its output); it writes the inputs of differential expression'.format(LEVELS_TABLE, Path.cwd()))
LEVELS = [line.split('\t') for line in LEVELS_TABLE.read_text().splitlines()[1:] if line]  # [column, level]
WHOLE = Path(DE_INPUT, 'model.tsv').is_file()  # absent with null comparisons only (#449); evaluate --de checks every set it declared
if not WHOLE and not LEVELS:
    raise WorkflowError('{} has no DE set: neither {} nor a level in {}; run `benchmark.py evaluate` in the project '
                        'first'.format(Path.cwd(), DE_INPUT, LEVELS_TABLE))
COMMAND = ('Rscript {params.script:q} --counts {input.counts:q} --library-sizes {input.library_sizes:q} --design {input.design:q} '
           '--model {input.model:q} --results {output.results:q} --record {output.record:q} > {log:q} 2>&1')


rule all:  # both outputs of every fit: Snakemake does not rerun a job for a missing output nobody requests
    input:
        [DE + '/' + name for name in ('results.tsv', 'record.tsv')] if WHOLE else [],
        ['{}/{}/{}/{}'.format(DE, column, level, name) for column, level in LEVELS for name in ('results.tsv', 'record.tsv')]


rule benchmark_de:
    input:
        counts = DE_INPUT + '/counts.tsv',
        library_sizes = DE_INPUT + '/library_sizes.tsv',
        design = DE_INPUT + '/design.tsv',
        model = DE_INPUT + '/model.tsv',
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
        model = DE + '/{column}/{level}/input/model.tsv',
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
