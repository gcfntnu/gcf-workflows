#-*- mode:snakemake -*-
"""Differential expression on a ratiometric pool's evaluation (#414): de.R in the `sartools` image of docker.config.

Run in the evaluated project after `benchmark.py evaluate` has written benchmark/de/input/, with the same container
settings as the project's `quant_all` (the repository holding de.R must be visible in the container):

    snakemake -s <repo>/smallrna/benchmark/de.smk -d <project> --use-apptainer --cores 1

then score the results with `benchmark.py evaluate --project <project> --de`.
"""
from pathlib import Path

import yaml

HERE = Path(str(workflow.current_basedir))
SARTOOLS = yaml.safe_load((HERE.parents[1] / 'docker.config').read_text())['docker']['sartools']
DE_INPUT = 'benchmark/de/input'


rule benchmark_de:
    input:
        counts = DE_INPUT + '/counts.tsv',
        library_sizes = DE_INPUT + '/library_sizes.tsv',
        design = DE_INPUT + '/design.tsv'
    output:
        results = 'benchmark/de/results.tsv',
        record = 'benchmark/de/record.tsv'
    params:
        script = str(HERE / 'de.R')
    threads:
        1
    resources:
        mem_mb = 2000
    log:
        'benchmark/de/de.log'
    container:
        'docker://' + SARTOOLS
    shell:
        'Rscript {params.script:q} --counts {input.counts:q} --library-sizes {input.library_sizes:q} --design {input.design:q} '
        '--results {output.results:q} --record {output.record:q} > {log:q} 2>&1'
