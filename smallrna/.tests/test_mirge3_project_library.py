"""mirge3_project_library.py: each project's miRge3 library links the shared one and has its own spike-in index.

bowtie-build is mocked.
"""
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/mirge3_project_library.py'
SUFFIXES = ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt')
CALIBRATORS = '>spikeIn|Cal01 spikeIn 01\nGTCCCACTCCGTAGATCTGTTC\n>spikeIn|Cal02 spikeIn 2\nGATGTAACGAGTTGGAATGCAA\n'
ERCC = '>ERCC-00002\nTCCAGATTACTTCCATTTCCGCCCAAGCTGCTCACAGTATACGGGCGTCGGCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTACG\n'

FAKE_BOWTIE_BUILD = '''#!/usr/bin/env python3
import shutil, sys
fasta, prefix = sys.argv[-2], sys.argv[-1]
for suffix in ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt'):
    shutil.copy(fasta, prefix + '.' + suffix)
'''


def snapshot(path):
    """Every file under path with its content, links not followed."""
    return {str(p.relative_to(path)): (os.readlink(p) if p.is_symlink() else p.read_text())
            for p in sorted(path.rglob('*')) if p.is_symlink() or p.is_file()}


@pytest.fixture
def setup(tmp_path):
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    (bin_dir / 'bowtie-build').write_text(FAKE_BOWTIE_BUILD)
    (bin_dir / 'bowtie-build').chmod(0o755)
    shared = tmp_path / 'ext/mirge3/09-11-2026/human'
    for name in ['index.Libs/human_mirna_miRBase.' + s for s in SUFFIXES] + \
            ['annotation.Libs/human_miRBase.gff3', 'fasta.Libs/human_mature_miRBase.fa', 'gcf_manifest.tsv']:
        (shared / name).parent.mkdir(parents=True, exist_ok=True)
        (shared / name).write_text('library ' + name)
    fasta = tmp_path / 'ext/spikein/fasta'
    fasta.mkdir(parents=True)
    (fasta / 'smallrna_calibrators.fa').write_text(CALIBRATORS)
    (fasta / 'ERCC92.fa').write_text(ERCC)
    env = dict(os.environ, PATH='{}:{}'.format(bin_dir, os.environ['PATH']))
    return tmp_path, shared, fasta, env


def run(setup, sets, project='p1'):
    tmp_path, shared, fasta, env = setup
    output = tmp_path / project / 'quant/mirge3/lib/human'
    record = tmp_path / project / 'quant/mirge3/spikein.tsv'
    result = subprocess.run([sys.executable, str(SCRIPT), '--library', str(shared), '--output', str(output),
                             '--record', str(record), *[str(fasta / (name + '.fa')) for name in sets]],
                            capture_output=True, text=True, env=env, timeout=60)
    rows = {row['field']: row['value'] for row in csv.DictReader(record.open(), delimiter='\t')} if record.exists() else {}
    return result, output, rows


def spikein_files(library):
    return sorted(path.name for path in (library / 'index.Libs').glob('human_spike-in.*'))


def test_project_library_links_the_shared_library_and_holds_its_own_spikein_index(setup):
    _, shared, _, _ = setup
    result, library, rows = run(setup, ['smallrna_calibrators'])
    assert result.returncode == 0, result.stderr
    assert (library / 'fasta.Libs').is_symlink() and os.readlink(library / 'fasta.Libs') == str(shared / 'fasta.Libs')
    index = library / 'index.Libs'
    assert index.is_dir() and not index.is_symlink()
    assert (index / 'human_mirna_miRBase.1.ebwt').read_text() == 'library index.Libs/human_mirna_miRBase.1.ebwt'
    assert spikein_files(library) == sorted('human_spike-in.' + s for s in SUFFIXES)
    assert not any(path.is_symlink() for path in index.glob('human_spike-in.*'))
    assert (index / 'human_spike-in.1.ebwt').read_text() == CALIBRATORS
    assert rows['sets'] == 'smallrna_calibrators' and rows['sequences'] == '2'
    assert rows['fasta_sha256'] == hashlib.sha256(CALIBRATORS.encode()).hexdigest()
    assert rows['library'] == str(shared)


def test_projects_keep_their_own_spikeins_and_the_shared_library_is_unchanged(setup):
    """Returning to an earlier project after another changed the selection finds that project's index unchanged."""
    _, shared, _, _ = setup
    before = snapshot(shared)
    _, first, _ = run(setup, ['smallrna_calibrators'], project='p1')
    _, second, _ = run(setup, ['ERCC92'], project='p2')
    _, third, rows = run(setup, [], project='p3')
    assert (first / 'index.Libs/human_spike-in.1.ebwt').read_text() == CALIBRATORS
    assert (second / 'index.Libs/human_spike-in.1.ebwt').read_text() == ERCC
    assert spikein_files(third) == [] and rows['status'] == 'no spike-in sets: no spike-in index'
    assert snapshot(shared) == before


def test_spikein_index_shipped_in_the_shared_library_is_not_linked(setup):
    _, shared, _, _ = setup
    (shared / 'index.Libs/human_spike-in.1.ebwt').write_text('shipped')
    result, library, _ = run(setup, [])
    assert result.returncode == 0, result.stderr
    assert spikein_files(library) == []


def test_both_sets_are_combined_in_the_given_order(setup):
    result, library, rows = run(setup, ['smallrna_calibrators', 'ERCC92'])
    assert result.returncode == 0, result.stderr
    assert rows['sets'] == 'smallrna_calibrators,ERCC92' and rows['sequences'] == '3'
    assert (library / 'index.Libs/human_spike-in.1.ebwt').read_text() == CALIBRATORS + ERCC


def test_rebuild_replaces_the_project_library(setup):
    run(setup, ['smallrna_calibrators'])
    result, library, rows = run(setup, [])
    assert result.returncode == 0, result.stderr
    assert spikein_files(library) == [] and rows['sets'] == 'none'
    assert [p.name for p in library.parent.iterdir()] == ['human']


def test_failed_build_leaves_no_project_library(setup):
    tmp_path, _, _, _ = setup
    (tmp_path / 'bin/bowtie-build').write_text('#!/bin/sh\ntouch "$4.1.ebwt"\necho broken >&2\nexit 1\n')
    result, library, rows = run(setup, ['smallrna_calibrators'])
    assert result.returncode != 0
    assert 'bowtie-build failed' in result.stderr and 'broken' in result.stderr
    assert rows == {} and not library.exists() and list(library.parent.iterdir()) == []


def test_incomplete_spikein_index_fails(setup):
    tmp_path, _, _, _ = setup
    (tmp_path / 'bin/bowtie-build').write_text('#!/bin/sh\ntouch "$4.1.ebwt" "$4.2.ebwt"\n')
    result, library, _ = run(setup, ['smallrna_calibrators'])
    assert result.returncode != 0
    assert 'human_spike-in.3.ebwt' in result.stderr and not library.exists()


def test_output_must_be_named_after_the_organism(setup):
    tmp_path, shared, _, env = setup
    result = subprocess.run([sys.executable, str(SCRIPT), '--library', str(shared), '--output', str(tmp_path / 'lib/mouse'),
                             '--record', str(tmp_path / 'r.tsv')], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode != 0 and "organism name 'human'" in result.stderr
