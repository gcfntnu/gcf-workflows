"""Synthetic level-2 delivery DAGs; no quantifiers, containers or downloads are run."""
from pathlib import Path
import os
import re
import subprocess
import sys

import pytest
import yaml
from configmaker.libprep import LibprepConfig, LibprepConfigError


WORKFLOW = Path(__file__).resolve().parents[1] / 'singlecell.smk'
LIBPREP = WORKFLOW.parents[1] / 'libprep.config'
V3_KIT = '10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3'
MULTIOME_KIT = '10X Genomics Chromium Single Cell Multiome Gene Expression'
PARSE_KIT = 'Parse Biosciences Evercode WT v4'


@pytest.mark.parametrize('kit,chemistry', [(V3_KIT, 'SC3Pv3'), (MULTIOME_KIT, 'ARC-v1')])
@pytest.mark.parametrize('read_geometry', [[28, 90], [28, 150]])
@pytest.mark.parametrize('suffix', ['', ' PE'])
def test_cellranger_libprep_selection(kit, chemistry, read_geometry, suffix):
    selection = LibprepConfig.load(LIBPREP).select(kit + suffix, read_geometry)
    assert selection.entry == kit + ' PE'
    assert selection.workflow == 'singlecell'
    assert selection.read_geometry == tuple(read_geometry)
    assert selection.parameters['quant']['method'] == 'cellranger'
    assert selection.parameters['quant']['cellranger']['chemistry'] == chemistry
    assert selection.parameters['db']['reference_db'] == '10xgenomics'
    assert selection.parameters['delta_readlen'] == [0, 0]
    if kit == MULTIOME_KIT:
        assert 'whitelist' not in selection.parameters['quant']['cellranger']
        assert not {'starsolo', 'alevin', 'umi_tools'} & selection.parameters['quant'].keys()
    else:
        assert selection.parameters['quant']['cellranger']['whitelist'] == '3M-february-2018_TRU.txt'


@pytest.mark.parametrize('kit', [V3_KIT, MULTIOME_KIT])
def test_cellranger_libprep_requires_paired_end(kit):
    with pytest.raises(LibprepConfigError, match='Unsupported library kit'):
        LibprepConfig.load(LIBPREP).select(kit, [28])


@pytest.mark.parametrize('method,kit,chemistry', [
    pytest.param('cellranger', V3_KIT, 'SC3Pv3', id='cellranger-v3'),
    pytest.param('cellranger', MULTIOME_KIT, 'ARC-v1', id='cellranger-multiome'),
    pytest.param('10x_starsolo', V3_KIT, None, id='10x-starsolo'),
    pytest.param('splitpipe', PARSE_KIT, None, id='splitpipe'),
    pytest.param('parsebio_starsolo', PARSE_KIT, None, id='parsebio-starsolo'),
])
@pytest.mark.parametrize('preprocess', [False, True])
def test_level2_delivery_dag(tmp_path, method, kit, chemistry, preprocess):
    parse = method in {'splitpipe', 'parsebio_starsolo'}
    config = {
        'workflow': 'singlecell', 'organism': 'mus_musculus', 'project_id': ['SYNTHETIC'],
        'libprepkit': kit,
        'machine': 'Synthetic sequencer', 'read_geometry': [100, 100] if parse else [28, 90],
        'ext_dir': 'data/ext', 'interim_dir': 'data/tmp', 'fastq_dir': 'data/raw/fastq',
        'samples': {sample: {'R1': f'{sample}_R1.fastq.gz', 'R2': f'{sample}_R2.fastq.gz'}
                    for sample in ['libA', 'libB']},
        'quant': {},
        'preprocessing': {'enabled': preprocess},
    }
    # Cell Ranger method, chemistry and aggregation must come from the preset/defaults.
    if method != 'cellranger':
        config['quant'].update(method=method, use_velo=False, aggregate={'method': 'scanpy'})
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
    if method == 'cellranger':
        assertions += f'''
assert config['quant']['method'] == BFQ_METHOD == AGGR_METHOD == 'cellranger'
assert SINGLECELL_ASSAY == 'wt'
assert config['quant']['cellranger']['chemistry'] == {chemistry!r}
assert config['read_geometry'] == [28, 90]
assert METHODS == ['cellranger']
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
    (tmp_path / 'dry-run.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    if method == 'cellranger':
        scheduled = re.findall(r'^(?:local)?rule (\w+):$', result.stdout, re.MULTILINE)
        assert not any(other in name for name in scheduled for other in ['starsolo', 'splitpipe', 'tcr', 'atac'])
        # Inspect each expanded count job: dependency tracking must retain both raw reads.
        jobs = re.findall(r'rule cellranger_quant_:\n(.*?cellranger count [^\n]*)', result.stdout, re.DOTALL)
        assert len(jobs) == 2
        for sample in config['samples']:
            job = next(job for job in jobs if f'wildcards: sample={sample}\n' in job)
            assert f'data/raw/fastq/{sample}_R1.fastq.gz' in job
            assert f'data/raw/fastq/{sample}_R2.fastq.gz' in job
            command = next(line for line in job.splitlines() if 'cellranger count ' in line)
            assert f'--chemistry {chemistry} ' in command
            assert '--fastqs data/raw/fastq ' in command
            assert f'--sample {sample} ' in command
        assert 'rule cellranger_aggr:' in result.stdout
        assert 'rule singlecell_metadata_tables:' in result.stdout
    if preprocess:
        assert 'rule preprocess_final_umap_png:' in result.stdout
    if method in {'10x_starsolo', 'parsebio_starsolo'}:
        assert 'rule bfq_level2_starsolo_aggr_mtx:' in result.stdout
