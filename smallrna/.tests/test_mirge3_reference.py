"""mirge3_reference.py: the selected miRge3 library holds the requested miRNA database, and what it records."""
import csv
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/mirge3_reference.py'
FILES = ['annotation.Libs/mouse_{db}.gff3', 'annotation.Libs/mouse_merges_{db}.csv',
         'annotation.Libs/mouse_miRNAs_in_repetitive_element_{db}.csv', 'fasta.Libs/mouse_mature_{db}.fa',
         'fasta.Libs/mouse_mirna_SNP_pseudo_{db}.fa', 'index.Libs/mouse_hairpin_{db}.1.ebwt', 'index.Libs/mouse_mirna_{db}.1.ebwt']


def make_library(tmp_path, databases=('miRBase', 'MirGeneDB'), provenance=True):
    lib = tmp_path / 'ext/mirge3/09-11-2026/mouse'
    for db in databases:
        for name in FILES:
            path = lib / name.format(db=db)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('synthetic')
    if provenance:
        (lib / 'gcf_mirna_db.tsv').write_text('database\trelease\nmirbase\t23.1\nmirgenedb\t3.0\n')
        (lib / 'gcf_manifest.tsv').write_text('key\tvalue\nensembl_release\t116\ngenome_assembly\tGRCm39\n')
        (lib / 'gcf_sources.tsv').write_text('file\turl\tsha256\tretrieved\na.fa\thttps://example.org/a.fa\t{}\tNA\n'.format('a' * 64))
    return lib


def run(tmp_path, lib, mirna_db='mirbase', verification=True):
    args = [sys.executable, str(SCRIPT), '--library', str(lib), '--mirna-db', mirna_db, '--release', '09-11-2026']
    if verification:
        record = tmp_path / 'verify.tsv'
        record.write_text('field\tvalue\nstatus\tverified\npublished_sha256\t{}\n'.format('b' * 64))
        args += ['--verification', str(record)]
    output = tmp_path / 'reference.tsv'
    result = subprocess.run(args + ['--output', str(output)], capture_output=True, text=True, timeout=60)
    rows = {row['field']: (row['value'], row['source']) for row in csv.DictReader(output.open(), delimiter='\t')} if output.exists() else {}
    return result, rows


@pytest.mark.parametrize('mirna_db, release', [('mirbase', '23.1'), ('mirgenedb', '3.0')])
def test_requested_database_and_its_release_are_recorded(tmp_path, mirna_db, release):
    result, rows = run(tmp_path, make_library(tmp_path), mirna_db)
    assert result.returncode == 0, result.stderr
    assert rows['mirna_database'] == (mirna_db, 'quant.mirge3.db')
    assert rows['mirna_release'] == (release, 'gcf_mirna_db.tsv')
    assert rows['ensembl_release'] == ('116', 'gcf_manifest.tsv')
    assert rows['source_checksums'] == ('1 files with SHA-256', 'gcf_sources.tsv')
    assert rows['archive_status'] == ('verified', 'verify.tsv')


def test_library_without_the_requested_database_fails_naming_the_files(tmp_path):
    result, rows = run(tmp_path, make_library(tmp_path, databases=('miRBase',)), 'mirgenedb')
    assert result.returncode != 0
    assert 'quant.mirge3.db is mirgenedb' in result.stderr
    assert 'index.Libs/mouse_mirna_MirGeneDB.1.ebwt' in result.stderr
    assert rows == {}


def test_provenance_not_listing_the_database_fails(tmp_path):
    lib = make_library(tmp_path)
    (lib / 'gcf_mirna_db.tsv').write_text('database\trelease\nmirbase\t23.1\n')
    result, _ = run(tmp_path, lib, 'mirgenedb')
    assert result.returncode != 0
    assert 'lists mirbase, not quant.mirge3.db mirgenedb' in result.stderr


def test_library_without_records_is_reported_as_unknown(tmp_path):
    """The SourceForge library has no provenance files and no verification record."""
    result, rows = run(tmp_path, make_library(tmp_path, provenance=False), verification=False)
    assert result.returncode == 0, result.stderr
    assert rows['mirna_release'] == ('NA', 'NA')
    assert rows['ensembl_release'] == ('NA', 'NA')
    assert rows['archive_status'][0] == 'not verified: library without published checksums'
