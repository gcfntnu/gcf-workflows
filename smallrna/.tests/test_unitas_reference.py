"""unitas_reference.py: what the selected UNITAS dump records, and its compatibility with the miRNA annotation."""
import csv
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/unitas_reference.py'
SERVED = '20,21,22,23.1'


def versions(mirna):
    return ('Genomic tRNA database version/date: ....... 05.11.2026 (dd.mm.yyyy)\n'
            'Ensembl version/date: ..................... Release 116\n'
            'tRF-1 sequence data version/date: ......... 09.04.2019 (dd.mm.yyyy)\n'
            'SILVA rRNA (SSU) database version/date: ... Release 132\n'
            'miRBase database version/date: ............ {}\n'.format(mirna))


def gff(database, release):
    return '##gff-version 3\n# microRNAs:               {} v{}\n# genome-build-id:         GRCm39\n'.format(database, release)


def make_dump(tmp_path, mirna='Release 23.1', mirna_tsv=None, sources=True):
    dump = tmp_path / 'ext/unitas/05-11-2026/UNITAS_refdump_mus_musculus'
    dump.mkdir(parents=True)
    if mirna is not None:
        (dump / 'db_versions.info').write_text(versions(mirna))
    if mirna_tsv:
        (dump / 'gcf_mirna_db.tsv').write_text('database\trelease\n{}\t{}\n'.format(*mirna_tsv))
    if sources:
        (dump / 'gcf_sources.tsv').write_text('file\turl\tsha256\tretrieved\nmature.fa\thttps://example.org/mature.fa\t{}\tNA\n'.format('a' * 64))
    return dump


def run(tmp_path, dump, annotation_text, mirna_db='mirbase', status='verified'):
    verification = tmp_path / 'verify.tsv'
    verification.write_text('field\tvalue\nstatus\t{}\npublished_sha256\t{}\nchecked\t2026-11-11T00:00:00Z\n'.format(
        status, 'b' * 64 if status == 'verified' else 'NA'))
    annotation = tmp_path / 'genes.gff3'
    annotation.write_text(annotation_text)
    output = tmp_path / 'reference.tsv'
    result = subprocess.run([sys.executable, str(SCRIPT), '--dump', str(dump), '--verification', str(verification),
                             '--mirna-db', mirna_db, '--annotation', str(annotation), '--release', '05-11-2026',
                             '--mirbase-releases', SERVED, '--output', str(output)], capture_output=True, text=True, timeout=60)
    record = {row['field']: (row['value'], row['source']) for row in csv.DictReader(output.open(), delimiter='\t')} if output.exists() else {}
    return result, record


@pytest.mark.parametrize('line', ['Release 23.1', '23.1'])
def test_matching_mirbase_release_is_recorded_from_unitas_versions(tmp_path, line):
    result, record = run(tmp_path, make_dump(tmp_path, mirna=line), gff('miRBase', '23.1'))
    assert result.returncode == 0, result.stderr
    assert record['mirna_database'] == ('mirbase', 'db_versions.info')
    assert record['mirna_release'] == ('23.1', 'db_versions.info')
    assert record['annotation_release'] == ('23.1', 'genes.gff3')
    assert record['ensembl_release'] == ('116', 'db_versions.info')
    assert record['gtrnadb'] == ('2026-11-05', 'db_versions.info')
    assert record['silva_ssu'] == ('132', 'db_versions.info')
    assert record['archive_status'] == ('verified', 'verify.tsv')
    assert record['source_checksums'] == ('1 files with SHA-256', 'gcf_sources.tsv')


def test_missing_records_are_reported_as_unknown_not_filled_in(tmp_path):
    result, record = run(tmp_path, make_dump(tmp_path, sources=False), gff('miRBase', '23.1'),
                         status='not verified: the release publishes no SHA256SUMS')
    assert result.returncode == 0, result.stderr
    assert record['source_checksums'] == ('NA', 'NA')
    assert record['seqmap'] == ('NA', 'NA')
    assert record['archive_status'][0] == 'not verified: the release publishes no SHA256SUMS'
    assert record['archive_published_sha256'][0] == 'NA'


def test_mirbase_release_mismatch_fails_naming_both_and_the_settings(tmp_path):
    result, record = run(tmp_path, make_dump(tmp_path, mirna='Release 22'), gff('miRBase', '23.1'))
    assert result.returncode != 0
    assert 'was built from miRBase 22, but the miRNA feature annotation' in result.stderr
    assert 'is miRBase 23.1' in result.stderr
    assert 'select a dump built from miRBase 23.1 (db.unitas.release) or set db.mirbase.release to 22' in result.stderr
    assert record == {}


def test_mismatch_with_a_release_gcfdb_cannot_annotate_does_not_suggest_it(tmp_path):
    result, _ = run(tmp_path, make_dump(tmp_path, mirna='22.1'), gff('miRBase', '23.1'))
    assert result.returncode != 0
    assert 'no annotation for miRBase 22.1 is available (db.mirbase.release: 20, 21, 22, 23.1)' in result.stderr
    assert 'set db.mirbase.release to 22.1' not in result.stderr


def test_mirgenedb_dump_selected_as_mirbase_fails(tmp_path):
    dump = make_dump(tmp_path, mirna='MirGeneDB 3.0', mirna_tsv=('mirgenedb', '3.0'))
    result, _ = run(tmp_path, dump, gff('miRBase', '23.1'))
    assert result.returncode != 0
    assert 'db.unitas.mirna_db is mirbase, but the selected dump' in result.stderr
    assert 'was built from MirGeneDB 3.0' in result.stderr


def test_mirgenedb_dump_with_its_own_gff(tmp_path):
    dump = make_dump(tmp_path, mirna='MirGeneDB 3.0', mirna_tsv=('mirgenedb', '3.0'))
    result, record = run(tmp_path, dump, gff('MirGeneDB', '3.0'), mirna_db='mirgenedb')
    assert result.returncode == 0, result.stderr
    assert record['mirna_database'] == ('mirgenedb', 'gcf_mirna_db.tsv')
    assert record['annotation_database'] == ('mirgenedb', 'genes.gff3')


def test_mirgenedb_dump_and_its_gff_disagreeing_fails(tmp_path):
    dump = make_dump(tmp_path, mirna='MirGeneDB 3.0', mirna_tsv=('mirgenedb', '3.0'))
    result, _ = run(tmp_path, dump, gff('MirGeneDB', '2.1'), mirna_db='mirgenedb')
    assert result.returncode != 0
    assert 'the dump and its gcf_mirgenedb.gff3 disagree' in result.stderr


@pytest.mark.parametrize('mirna, mirna_tsv, message', [
    ('Release 23.1', ('mirgenedb', '3.0'), 'records MirGeneDB 3.0, but db_versions.info records miRBase 23.1'),
    (None, None, 'cannot read the miRNA database release'),
    ('unknown', None, 'cannot read the miRNA database release'),
])
def test_inconsistent_or_unreadable_dump_records_fail(tmp_path, mirna, mirna_tsv, message):
    result, _ = run(tmp_path, make_dump(tmp_path, mirna=mirna, mirna_tsv=mirna_tsv), gff('miRBase', '23.1'))
    assert result.returncode != 0
    assert message in result.stderr


def test_annotation_without_release_header_fails(tmp_path):
    result, _ = run(tmp_path, make_dump(tmp_path), '##gff-version 3\n')
    assert result.returncode != 0
    assert 'has no "# microRNAs: miRBase v<release>"' in result.stderr
