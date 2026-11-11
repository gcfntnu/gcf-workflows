"""mirge3_spikein.py: the miRge3 spike-in index follows each project's spike-in sets; bowtie-build is mocked."""
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/mirge3_spikein.py'
CALIBRATORS = '>spikeIn|Cal01 spikeIn 01\nGTCCCACTCCGTAGATCTGTTC\n>spikeIn|Cal02 spikeIn 2\nGATGTAACGAGTTGGAATGCAA\n'
ERCC = '>ERCC-00002\nTCCAGATTACTTCCATTTCCGCCCAAGCTGCTCACAGTATACGGGCGTCGGCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTACGAGCACTATGGTCAGTAATTCCTGGAGGAATAGGTACCAAGAAAAAAACGAACCTTTGGGTTCCAGAGCTGTACGGTCGCACTGAACTCGGATAGGTCTCAGAAAAACGAAATATAGGCTTACGGTAGGTCCGAATGGCACAAAGCTTGTTCCGTTAGCTGGCATAAGATTCCATGCCTAGATGTGGTACACGTTTTCGCAGCCGTGCCATGTGTAGATCATTGAACCAGTGTGGGCAA\n'

FAKE_BOWTIE_BUILD = '''#!/usr/bin/env python3
import shutil, sys
fasta, prefix = sys.argv[-2], sys.argv[-1]
for suffix in ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt'):
    shutil.copy(fasta, prefix + '.' + suffix)
'''


@pytest.fixture
def setup(tmp_path):
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    (bin_dir / 'bowtie-build').write_text(FAKE_BOWTIE_BUILD)
    (bin_dir / 'bowtie-build').chmod(0o755)
    index = tmp_path / 'ext/mirge3/09-11-2026/human/index.Libs'
    index.mkdir(parents=True)
    (index / 'human_mirna_miRBase.1.ebwt').write_text('library index')
    fasta = tmp_path / 'ext/spikein/fasta'
    fasta.mkdir(parents=True)
    (fasta / 'smallrna_calibrators.fa').write_text(CALIBRATORS)
    (fasta / 'ERCC92.fa').write_text(ERCC)
    env = dict(os.environ, PATH='{}:{}'.format(bin_dir, os.environ['PATH']))
    return tmp_path, index, fasta, env


def run(setup, sets, project='p1'):
    tmp_path, index, fasta, env = setup
    output = tmp_path / project / 'quant/mirge3/spikein.tsv'
    result = subprocess.run([sys.executable, str(SCRIPT), '--index-dir', str(index), '--organism-name', 'human',
                             '--output', str(output), *[str(fasta / (name + '.fa')) for name in sets]],
                            capture_output=True, text=True, env=env, timeout=60)
    record = {row['field']: row['value'] for row in csv.DictReader(output.open(), delimiter='\t')} if output.exists() else {}
    return result, record


def index_files(index):
    return sorted(path.name for path in index.glob('human_spike-in.*'))


def test_index_is_built_from_the_configured_sets_only(setup):
    _, index, _, _ = setup
    result, record = run(setup, ['smallrna_calibrators'])
    assert result.returncode == 0, result.stderr
    assert record['sets'] == 'smallrna_calibrators'
    assert record['sequences'] == '2'
    assert record['fasta_sha256'] == hashlib.sha256(CALIBRATORS.encode()).hexdigest()
    assert len(index_files(index)) == 6
    assert (index / 'human_spike-in.1.ebwt').read_text() == CALIBRATORS
    assert (index / 'human_mirna_miRBase.1.ebwt').read_text() == 'library index'


def test_next_project_replaces_the_index_with_its_own_sets(setup):
    _, index, _, _ = setup
    run(setup, ['smallrna_calibrators'], project='p1')
    (index / 'human_spike-in.stale.ebwt').write_text('left from an older build')
    result, record = run(setup, ['ERCC92'], project='p2')
    assert result.returncode == 0, result.stderr
    assert record['sets'] == 'ERCC92'
    assert (index / 'human_spike-in.1.ebwt').read_text() == ERCC
    assert 'human_spike-in.stale.ebwt' not in index_files(index)


def test_both_sets_are_combined_in_the_given_order(setup):
    _, index, _, _ = setup
    result, record = run(setup, ['smallrna_calibrators', 'ERCC92'])
    assert result.returncode == 0, result.stderr
    assert record['sets'] == 'smallrna_calibrators,ERCC92'
    assert record['sequences'] == '3'
    assert (index / 'human_spike-in.1.ebwt').read_text() == CALIBRATORS + ERCC


def test_no_sets_removes_the_index(setup):
    _, index, _, _ = setup
    run(setup, ['ERCC92'], project='p1')
    result, record = run(setup, [], project='p2')
    assert result.returncode == 0, result.stderr
    assert record['sets'] == 'none'
    assert record['status'] == 'no spike-in sets: index removed'
    assert index_files(index) == []
    assert (index / 'human_mirna_miRBase.1.ebwt').exists()


def test_failed_build_leaves_no_partial_index(setup):
    tmp_path, index, _, env = setup
    (tmp_path / 'bin/bowtie-build').write_text('#!/bin/sh\ntouch "$4.1.ebwt"\necho broken >&2\nexit 1\n')
    result, record = run(setup, ['smallrna_calibrators'])
    assert result.returncode != 0
    assert 'bowtie-build failed' in result.stderr and 'broken' in result.stderr
    assert record == {}
    assert index_files(index) == []
