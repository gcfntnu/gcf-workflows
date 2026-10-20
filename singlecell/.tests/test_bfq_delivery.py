"""Synthetic level-2 delivery DAGs; no quantifiers, containers or downloads are run."""
from pathlib import Path
import os
import subprocess
import sys

import pytest
import yaml


WORKFLOW = Path(__file__).resolve().parents[1] / 'singlecell.smk'


@pytest.mark.parametrize('method', ['cellranger', '10x_starsolo', 'splitpipe', 'parsebio_starsolo'])
@pytest.mark.parametrize('preprocess', [False, True])
def test_level2_delivery_dag(tmp_path, method, preprocess):
    parse = method in {'splitpipe', 'parsebio_starsolo'}
    config = {
        'workflow': 'singlecell', 'organism': 'mus_musculus', 'project_id': ['SYNTHETIC'],
        'libprepkit': 'Parse Biosciences Evercode WT v4' if parse else
            '10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3',
        'machine': 'Synthetic sequencer', 'read_geometry': [100, 100],
        'ext_dir': 'data/ext', 'interim_dir': 'data/tmp', 'fastq_dir': 'data/raw/fastq',
        'samples': {sample: {'R1': f'{sample}_R1.fastq.gz', 'R2': f'{sample}_R2.fastq.gz'}
                    for sample in ['libA', 'libB']},
        'quant': {'method': method, 'use_velo': False,
                  'aggregate': {'method': 'cellranger' if method == 'cellranger' else 'scanpy'}},
        'preprocessing': {'enabled': preprocess},
    }
    if parse:
        config['wells'] = {'bioA': {'Wells': 'A1'}, 'bioB': {'Wells': 'B1'}}
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(config))
    (tmp_path / 'pep.yaml').write_text('pep_version: 2.0.0\nsample_table: samples.csv\n')
    (tmp_path / 'samples.csv').write_text('sample_name,R1,R2,paired_end\n'
        'libA,libA_R1.fastq.gz,libA_R2.fastq.gz,1\nlibB,libB_R1.fastq.gz,libB_R2.fastq.gz,1\n')
    fastqs = tmp_path / 'data/raw/fastq'
    fastqs.mkdir(parents=True)
    for sample in config['samples']:
        for read in ['R1', 'R2']:
            (fastqs / f'{sample}_{read}.fastq.gz').touch()

    # Evaluate the actual workflow objects and then construct its entire reporting DAG.
    assertions = '''
assert 'data/tmp/singlecell/bfq/exprs/scanpy/all_samples_filtered.h5ad' in rules.bfq_level2_exprs.output
assert ('data/tmp/singlecell/bfq/exprs/scanpy/all_samples_preprocessed.h5ad' in rules.bfq_level2_exprs.output) == PREPROCESS_ENABLED
assert bool(BFQ_PREPROCESS_FIGS) == PREPROCESS_ENABLED
if PREPROCESS_ENABLED:
    assert Path(rules.preprocess_final_umap_png.params.script).is_file()
if BFQ_METHOD == 'cellranger':
    produced = set(expand(rules.cellranger_aggr.output, aggr_id=AGGR_IDS))
    assert set(rules.bfq_level2_logs.input.aggr_summary).issubset(produced)
    assert set(rules.bfq_level2_logs.input.aggr_csv).issubset(produced)
if BFQ_METHOD in {'10x_starsolo', 'parsebio_starsolo'}:
    assert Path(rules.bfq_level2_starsolo_aggr_mtx.params.script).is_file()
    for name in ['matrix.mtx', 'features.tsv', 'barcodes.tsv']:
        assert f'data/tmp/singlecell/bfq/exprs/mtx/{BFQ_METHOD}/all_samples/{name}' in rules.bfq_all.input
'''
    (tmp_path / 'Snakefile').write_text('from pathlib import Path\npepfile: "pep.yaml"\nconfigfile: "config.yaml"\n'
        + f'include: {str(WORKFLOW)!r}\n' + assertions)
    env = os.environ.copy()
    env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env['PATH']
    for name in ['GCF_SECRET', 'GCF_EXT', 'GCF_INTERIM', 'GCF_FASTQ', 'GCF_REF']:
        env.pop(name, None)
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '-p', '--cores', '1',
        'multiqc_report'], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    if preprocess:
        assert 'rule preprocess_final_umap_png:' in result.stdout
    if method in {'10x_starsolo', 'parsebio_starsolo'}:
        assert 'rule bfq_level2_starsolo_aggr_mtx:' in result.stdout
