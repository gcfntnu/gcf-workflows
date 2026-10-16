"""Synthetic Parse metadata contracts; no scientific tools or facility inputs."""
import copy
import importlib.util
from pathlib import Path
import subprocess
import sys

import anndata
import numpy as np
import pandas as pd
import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[1] / 'rules/quant/scripts'


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


metadata = load('create_metadata_tables')
convert = load('convert_scanpy')


def config():
    return {'libprepkit': 'Parse Biosciences Test', 'samples': {
        '001': {'Sample_ID': '001', 'Concentration': 0, 'custom': 'lib A', 'R1': 'a.fastq.gz', 'I1': 'a-index.fastq.gz'},
        '1': {'Sample_ID': '1', 'Concentration': 3, 'custom': 'lib B', 'R1': 'b.fastq.gz', 'I1': 'b-index.fastq.gz'},
    }, 'wells': {
        '007': {'Sample_ID': '007', 'Wells': 'A1-A2', 'custom': 'bio A', 'External_ID': '0007', 'zero': 0},
        '7': {'Sample_ID': '7', 'Wells': 'B1', 'custom': 'bio B', 'External_ID': '0008', 'zero': 0},
    }}


@pytest.mark.parametrize('equal', [False, True])
def test_section_ownership_and_values(equal):
    source = config()
    if equal:
        source['samples']['1']['Concentration'] = 0
    before = copy.deepcopy(source)
    samples, libraries = metadata.normalized_tables(source)
    assert source == before
    assert samples.Sample_ID.tolist() == ['007', '7']
    assert libraries.library_id.tolist() == ['001', '1']
    assert samples.Wells.tolist() == ['A1-A2', 'B1']
    assert samples.zero.tolist() == [0, 0]
    assert libraries.Concentration.tolist() == [0, 0 if equal else 3]
    assert 'Concentration' not in samples and 'Sample_ID' not in libraries
    assert libraries.I1.tolist() == ['a-index.fastq.gz', 'b-index.fastq.gz']


@pytest.mark.parametrize('field', ['Sample_Group', 'Sample_Biosource', 'External_ID'])
@pytest.mark.parametrize('value', ['control', 0, False])
def test_direct_config_rejects_biological_library_fields(field, value):
    source = config()
    source['samples']['001'][field] = value
    with pytest.raises(ValueError, match='Cell Multiplexing'):
        metadata.normalized_tables(source)


@pytest.mark.parametrize('value', [None, '', ' ', 'NA', float('nan')])
def test_empty_reserved_fields_are_not_library_annotations(value):
    source = config()
    source['samples']['001']['External_ID'] = value
    _, libraries = metadata.normalized_tables(source)
    assert 'External_ID' not in libraries


@pytest.mark.parametrize('section,field', [('samples', 'Sample_ID'), ('samples', 'library_id'), ('wells', 'Sample_ID')])
def test_conflicting_embedded_identity_fails(section, field):
    source = config()
    first = next(iter(source[section].values()))
    first[field] = 'WRONG'
    with pytest.raises(ValueError, match='disagrees|conflicting'):
        metadata.normalized_tables(source)


def test_exact_whitespace_identity_is_not_trimmed():
    source = config()
    source['samples'][' 001 '] = {'Sample_ID': ' 001 ', 'custom': 'value'}
    _, libraries = metadata.normalized_tables(source)
    assert libraries.library_id.tolist() == ['001', '1', ' 001 ']


def test_non_parse_split_unchanged():
    source = config()
    source['libprepkit'] = '10x Test'
    samples, libraries = metadata.normalized_tables(source)
    assert samples.Sample_ID.tolist() == ['001', '1']
    assert 'Concentration' in samples and 'Concentration' not in libraries
    assert 'R1' in libraries


def tables(tmp_path):
    source = tmp_path / 'config.yaml'
    source.write_text(yaml.safe_dump(config()))
    sample, library = tmp_path / 'sample_info.tsv', tmp_path / 'library_info.tsv'
    subprocess.run([sys.executable, str(SCRIPTS / 'create_metadata_tables.py'), '--configfile', str(source),
                    '--sample-info', str(sample), '--library-info', str(library)], check=True)
    barcode = tmp_path / 'barcode_info.tsv'
    barcode.write_text('barcode\tSample_ID\tlibrary_id\ncell1\t007\t001\ncell2\t7\t001\ncell3\t007\t1\ncell4\t7\t1\n')
    return sample, library, barcode


def test_tsv_parser_broadcast_and_h5ad_roundtrip(tmp_path):
    sample, library, barcode = tables(tmp_path)
    args = convert.create_parser().parse_args(['-f', 'splitpipe', '-o', str(tmp_path / 'out.h5ad'),
        '--sample-info', str(sample), '--library-info', str(library), '--barcode-info', str(barcode)])
    assert args.parse_metadata
    assert set(args.sample_info.index) == {'007', '7'}
    assert set(args.library_info.index) == {'001', '1'}
    obs = args.barcode_info[0]
    obs = convert.anndata_friendly_dtypes(obs, protect_cols=('Sample_ID', 'library_id'))
    obs = convert.broadcast_entity_metadata(obs, args.sample_info, 'Sample_ID', 'sample_info')
    lib = convert._drop_blacklisted_columns(args.library_info, convert._LIBRARY_INFO_BLACKLIST)
    lib = convert.parse_library_columns(lib, args.sample_info)
    obs = convert.broadcast_entity_metadata(obs, lib, 'library_id', 'library_info')
    assert obs.Sample_ID.tolist() == ['007', '7', '007', '7']
    assert obs.library_id.tolist() == ['001', '001', '1', '1']
    assert obs.custom.tolist() == ['bio A', 'bio B', 'bio A', 'bio B']
    assert obs.library_custom.tolist() == ['lib A', 'lib A', 'lib B', 'lib B']
    assert obs.Concentration.tolist() == ['0', '0', '3', '3']
    assert obs.External_ID.tolist() == ['0007', '0008', '0007', '0008']
    assert 'Wells' not in obs and 'R1' not in obs
    data = anndata.AnnData(X=np.ones((4, 2)), obs=obs)
    data.write_h5ad(tmp_path / 'out.h5ad')
    restored = anndata.read_h5ad(tmp_path / 'out.h5ad')
    assert restored.obs.Sample_ID.tolist() == obs.Sample_ID.tolist()
    assert restored.obs.library_custom.tolist() == obs.library_custom.tolist()


@pytest.mark.parametrize('extra', ['library_custom', 'LIBRARY_CUSTOM'])
def test_generated_prefix_collision_fails(extra):
    samples, libraries = metadata.normalized_tables(config())
    libraries[extra] = 'already used'
    with pytest.raises(ValueError, match='prefix collides'):
        convert.parse_library_columns(libraries.set_index('library_id'), samples.set_index('Sample_ID'))


def test_equal_cross_level_column_is_still_preserved_twice():
    samples = pd.DataFrame({'custom': ['same']}, index=['bio'])
    libraries = pd.DataFrame({'custom': ['same']}, index=['lib'])
    assert list(convert.parse_library_columns(libraries, samples)) == ['library_custom']


@pytest.mark.parametrize('failure', ['duplicate', 'missing'])
def test_broadcast_reference_checks(failure):
    obs = pd.DataFrame({'library_id': ['unknown']}, index=['cell'])
    libraries = pd.DataFrame({'value': [1, 2]}, index=['lib', 'lib' if failure == 'duplicate' else 'other'])
    with pytest.raises(ValueError, match='not unique|missing from metadata'):
        convert.broadcast_entity_metadata(obs, libraries, 'library_id', 'library_info')


@pytest.mark.parametrize('transpose', [False, True], ids=['cells-by-genes', 'genes-by-cells'])
def test_splitpipe_reader_supports_source_only_mmread(tmp_path, monkeypatch, transpose):
    import logging
    from scipy.io import mmread, mmwrite
    from scipy.sparse import csr_matrix, isspmatrix_csr

    sample, library, barcode = tables(tmp_path)
    directory = tmp_path / '007' / 'DGE_filtered'
    directory.mkdir(parents=True)
    matrix = directory / 'DGE.mtx'
    expected = np.array([[1, 0, 4], [0, 0, 7]], dtype=np.int64)
    mmwrite(matrix, csr_matrix(expected.T if transpose else expected))
    (directory / 'all_genes.csv').write_text('gene_id,gene_name\ng1,GeneA\ng2,\ng3,GeneC\n')
    # The canonical mapping has a different order and also contains other samples.
    (directory / 'cell_metadata.csv').write_text('bc_wells\ncell3\ncell1\n')
    args = convert.create_parser().parse_args(['-f', 'splitpipe', '-o', str(tmp_path / 'out.h5ad'),
        '--barcode-rename', 'skip', '--sample-info', str(sample), '--library-info', str(library),
        '--barcode-info', str(barcode)])

    # SciPy 1.11 accepts only source; use real Matrix Market I/O behind that API.
    def source_only_mmread(source):
        return mmread(source)

    monkeypatch.setattr(convert, 'mmread', source_only_mmread)
    monkeypatch.setattr(convert, 'logger', logging.getLogger(__name__), raising=False)
    monkeypatch.setattr(convert, '_USE_VELO', False)
    data = convert.read_splitpipe(str(matrix), args)

    assert isspmatrix_csr(data.X)
    assert np.issubdtype(data.X.dtype, np.integer)
    np.testing.assert_array_equal(data.X.toarray(), expected)
    assert data.obs_names.tolist() == ['cell3', 'cell1']
    assert data.var_names.tolist() == ['g1', 'g2', 'g3']
    assert data.var.gene.tolist() == ['GeneA', 'g2', 'GeneC']
    assert data.X[:, 1].nnz == 0  # Retain the zero-count feature and its metadata.


@pytest.mark.parametrize('input_format', ['splitpipe', 'splitpipe_aggr'])
def test_converter_cli_preserves_two_level_metadata(tmp_path, input_format):
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix

    sample, library, barcode = tables(tmp_path)
    matrices = []
    for sample_id, cells in [('007', ['cell1', 'cell3']), ('7', ['cell2', 'cell4'])]:
        directory = tmp_path / sample_id / 'DGE_filtered'
        directory.mkdir(parents=True)
        matrix = directory / 'DGE.mtx'
        mmwrite(matrix, csr_matrix([[1, 2, 0], [0, 3, 1]]))
        (directory / 'all_genes.csv').write_text('gene_id,gene_name\ng1,GeneA\ng2,GeneB\ng3,GeneC\n')
        (directory / 'cell_metadata.csv').write_text('bc_wells\n' + '\n'.join(cells) + '\n')
        matrices.append(str(matrix))
    output = tmp_path / 'filtered.h5ad'
    result = subprocess.run([sys.executable, str(SCRIPTS / 'convert_scanpy.py'), *matrices,
        '-f', input_format, '-o', str(output), '--canonical-filtered', '--barcode-rename', 'skip',
        '--sample-info', str(sample), '--library-info', str(library), '--barcode-info', str(barcode)],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    data = anndata.read_h5ad(output)
    obs = data.obs.loc[['cell1', 'cell2', 'cell3', 'cell4']]
    assert obs.Sample_ID.tolist() == ['007', '7', '007', '7']
    assert obs.library_id.tolist() == ['001', '001', '1', '1']
    assert obs.custom.tolist() == ['bio A', 'bio B', 'bio A', 'bio B']
    assert obs.library_custom.tolist() == ['lib A', 'lib A', 'lib B', 'lib B']
    assert obs.External_ID.tolist() == ['0007', '0008', '0007', '0008']
    assert obs.Concentration.tolist() == ['0', '0', '3', '3']
    assert data.shape == (4, 3) and data.X.sum() == 14


def test_export_does_not_coerce_integer_beside_missing_value(tmp_path):
    source = config()
    del source['samples']['1']['Concentration']
    _, library = metadata.normalized_tables(source)
    path = tmp_path / 'library.tsv'
    metadata.write_table(library, path)
    read = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    assert read.Concentration.tolist() == ['0', '']


@pytest.mark.parametrize('method', ['splitpipe', 'parsebio_starsolo'])
def test_parse_workflow_metadata_delivery(tmp_path, method):
    import os

    source = config()
    # Quantifier sample-name normalization is outside this metadata adapter's scope.
    source['wells'] = {'bio' + key: {**value, 'Sample_ID': 'bio' + key} for key, value in source['wells'].items()}
    source.update(libprepkit='Parse Biosciences Evercode WT v3', workflow='singlecell',
                  organism='homo_sapiens', read_geometry=[100, 100], project_id=['SYNTHETIC'],
                  machine='Synthetic sequencer', quant={'method': method})
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(source))
    workflow = SCRIPTS.parents[2] / 'singlecell.smk'
    (tmp_path / 'pep.yaml').write_text('pep_version: 2.0.0\nsample_table: samples.csv\n')
    (tmp_path / 'samples.csv').write_text('sample_name,R1,R2,paired_end\n001,a.fastq.gz,a2.fastq.gz,1\n1,b.fastq.gz,b2.fastq.gz,1\n')
    (tmp_path / 'Snakefile').write_text(
        'pepfile: "pep.yaml"\nconfigfile: "config.yaml"\n' + f'include: {str(workflow)!r}\n' +
        'assert set(rules.bfq_parse_metadata.output).issubset(set(rules.bfq_all.input))\n' +
        'assert str(rules.sample_info.output[0]) == "data/tmp/sample_info.tsv"\n')
    env = os.environ.copy()
    env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env['PATH']
    command = [str(Path(sys.executable).with_name('snakemake')), '--cores', '1', '--printshellcmds', 'bfq_parse_metadata']
    for dryrun in [True, False]:
        result = subprocess.run(command + (['--dry-run'] if dryrun else []), cwd=tmp_path,
                                env=env, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
    for name in ['sample_info.tsv', 'library_info.tsv']:
        delivered = tmp_path / 'data/tmp/singlecell/bfq/metadata' / name
        original = tmp_path / 'data/tmp/singlecell/quant/metadata' / name
        assert delivered.read_bytes() == original.read_bytes()
        assert len(pd.read_csv(delivered, sep='\t')) == 2


def test_generated_prefix_cannot_reuse_observation_column():
    samples = pd.DataFrame({'custom': ['same']}, index=['bio'])
    libraries = pd.DataFrame({'custom': ['same']}, index=['lib'])
    with pytest.raises(ValueError, match='prefix collides'):
        convert.parse_library_columns(libraries, samples, ['library_custom'])


def lightweight_parse_fixture(tmp_path, method):
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix

    source = {'workflow': 'singlecell', 'organism': 'homo_sapiens', 'project_id': ['SYNTHETIC'],
        'libprepkit': 'Parse Biosciences Evercode WT v4', 'read_geometry': [64, 58],
        'machine': 'Synthetic sequencer', 'ref_dir': 'reference',
        'interim_dir': 'data/tmp', 'ext_dir': 'data/ext', 'fastq_dir': 'data/raw/fastq',
        'samples': {'lib1': {'R1': 'a.fastq.gz', 'R2': 'a2.fastq.gz'},
                    'lib2': {'R1': 'b.fastq.gz', 'R2': 'b2.fastq.gz'}},
        'wells': {'bioA': {'Wells': 'A1'}, 'bioB': {'Wells': 'B1'}},
        'quant': {'method': method, 'use_velo': False, 'starsolo': {'feature_count': 'Gene'}}}
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(source))
    (tmp_path / 'pep.yaml').write_text('pep_version: 2.0.0\nsample_table: samples.csv\n')
    (tmp_path / 'samples.csv').write_text(
        'sample_name,R1,R2,paired_end\nlib1,a.fastq.gz,a2.fastq.gz,1\nlib2,b.fastq.gz,b2.fastq.gz,1\n')
    workflow = SCRIPTS.parents[2] / 'singlecell.smk'
    (tmp_path / 'Snakefile').write_text(
        'pepfile: "pep.yaml"\nconfigfile: "config.yaml"\n' + f'include: {str(workflow)!r}\n')
    reference = tmp_path / 'reference/anno/genes.tsv'
    reference.parent.mkdir(parents=True)
    reference.write_text('gene_id\tgene_symbols\ng1\tGeneA\ng2\tGeneB\ng3\tGeneC\n')
    base = tmp_path / 'data/tmp/singlecell/quant/aggregate' / method / 'bioA'
    expected = np.array([[1, 2, 0], [0, 3, 4]], dtype=np.int64)
    barcodes = ['02_01_01__s2', '01_01_01__s1']
    if method == 'splitpipe':
        matrix = base / 'DGE_filtered/count_matrix.mtx'
        matrix.parent.mkdir(parents=True)
        mmwrite(matrix, csr_matrix(expected))
        (matrix.parent / 'all_genes.csv').write_text('gene_id,gene_name\ng1,GeneA\ng2,GeneB\ng3,GeneC\n')
        (matrix.parent / 'cell_metadata.csv').write_text(
            'bc_wells,sample\n' + ''.join(barcode + ',bioA\n' for barcode in barcodes))
    else:
        matrix = base / 'Solo.out/Gene/filtered/matrix.mtx'
        matrix.parent.mkdir(parents=True)
        mmwrite(matrix, csr_matrix(expected.T))
        (matrix.parent / 'features.tsv').write_text('g1\tGeneA\tGene Expression\ng2\tGeneB\tGene Expression\ng3\tGeneC\tGene Expression\n')
        (matrix.parent / 'barcodes.tsv').write_text('\n'.join(barcodes) + '\n')
        (base / 'barcode_info.tsv').write_text(
            'barcode\tsource_barcode\tSample_ID\tlibrary_id\n'
            '02_01_01__s2\t02_01_01\tbioA\tlib2\n01_01_01__s1\t01_01_01\tbioA\tlib1\n')
    target = Path('_tmp') / method / 'filtered/bioA/anndata.light.h5ad'
    return target, base / 'barcode_info.tsv', expected, barcodes


@pytest.mark.parametrize('method', ['splitpipe', 'parsebio_starsolo', 'parsebio_starsolo_rt'])
def test_parse_lightweight_primary_metadata_delivery(tmp_path, method):
    import os
    import shlex
    from scipy.io import mmread
    from scipy.sparse import isspmatrix_csr

    target, primary, expected, barcodes = lightweight_parse_fixture(tmp_path, method)
    env = os.environ.copy()
    env.pop('GCF_SECRET', None)
    env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env['PATH']
    command = [str(Path(sys.executable).with_name('snakemake')), str(target), '--cores', '1',
               '--printshellcmds', '--allowed-rules', 'tmp_lightweight_filtered', 'splitpipe_barcode_info']
    dryrun = subprocess.run(command + ['--dry-run'], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert dryrun.returncode == 0, dryrun.stdout + dryrun.stderr
    output = dryrun.stdout + dryrun.stderr
    converter_line = next(line for line in output.splitlines() if 'convert_scanpy.py ' in line)
    argv = shlex.split(converter_line)
    primary_relative = str(primary.relative_to(tmp_path))
    assert argv[argv.index('--barcode-info') + 1] == primary_relative
    assert argv[argv.index('-f') + 1] == ('parsebio_starsolo' if method.startswith('parsebio') else method)
    assert any(line.strip().startswith('input:') and primary_relative in line for line in output.splitlines())
    assert 'auto_qc/' not in converter_line and 'doublet' not in converter_line
    if method != 'splitpipe':
        return  # Other Parse variants have DAG coverage; no STARsolo processing is run.

    # Generate the primary mapping with the real producer, then run the actual rule.
    assert not primary.exists()
    assert 'rule splitpipe_barcode_info:' in output
    result = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    mapping = pd.read_csv(primary, sep='\t', dtype=str)
    assert mapping.barcode.tolist() == barcodes
    assert mapping.Sample_ID.tolist() == ['bioA', 'bioA']
    assert mapping.library_id.tolist() == ['lib2', 'lib1']
    data = anndata.read_h5ad(tmp_path / target)
    assert isspmatrix_csr(data.X)
    np.testing.assert_array_equal(data.X.toarray(), expected)
    assert data.obs_names.tolist() == barcodes
    assert data.var_names.tolist() == ['g1', 'g2', 'g3']
    assert data.var.gene_symbol.tolist() == ['GeneA', 'GeneB', 'GeneC']
    exported = tmp_path / target.parent / 'anndata.mtx_v2'
    np.testing.assert_array_equal(mmread(exported / 'matrix.mtx').T.toarray(), expected)
    assert (exported / 'barcodes.tsv').read_text().splitlines() == barcodes


@pytest.mark.parametrize('damage', ['missing', 'wrong_barcode', 'wrong_sample'])
def test_splitpipe_lightweight_requires_exact_primary_mapping(tmp_path, damage):
    target, primary, _, _ = lightweight_parse_fixture(tmp_path, 'splitpipe')
    mapping = ('barcode\tsource_barcode\tSample_ID\tlibrary_id\n'
               '02_01_01__s2\t02_01_01\tbioA\tlib2\n01_01_01__s1\t01_01_01\tbioA\tlib1\n')
    if damage == 'wrong_barcode':
        mapping = mapping.replace('01_01_01__s1', '99_01_01__s1')
    elif damage == 'wrong_sample':
        mapping = mapping.replace('bioA', 'bioB')
    primary.write_text(mapping)
    command = [sys.executable, str(SCRIPTS / 'convert_scanpy.py'),
        str(primary.parent / 'DGE_filtered/count_matrix.mtx'), '-f', 'splitpipe',
        '--barcode-rename', 'skip', '-F', 'anndata_lightweight', 'v2_mtx',
        '-o', str(tmp_path / target.parent / 'anndata')]
    if damage != 'missing':
        command.extend(['--barcode-info', str(primary)])
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    expected_error = {'missing': 'requires primary barcode_info', 'wrong_barcode': 'canonical barcode mismatch',
                      'wrong_sample': 'no barcode_info rows found'}[damage]
    assert expected_error in result.stderr
    assert not (tmp_path / target).exists()
