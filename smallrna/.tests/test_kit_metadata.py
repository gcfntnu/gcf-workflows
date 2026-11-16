"""Small-RNA kit metadata (#361): read 1 orientation and the NEXTflex blockers recorded in the report header."""
from pathlib import Path

import pytest
import yaml

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
KITS = yaml.safe_load((ROOT / 'libprep.config').read_text())
SMALLRNA = {name: kit for name, kit in KITS.items() if kit['workflow'] == 'smallrna'}
V3 = 'Bioo Scientific NEXTflex Small RNA-Seq Kit v3'
V4 = 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4'


@pytest.mark.parametrize('name', SMALLRNA)
def test_read_1_is_sense(name):
    """Ligation libraries: the 5' adapter carries the read 1 primer site, so read 1 reads the RNA as it is."""
    assert SMALLRNA[name]['read_orientation'] == 'forward' and SMALLRNA[name]['library_strand'] == 'second strand'


def test_only_nextflex_v4_has_blockers_none_by_default():
    assert {name: kit['blockers'] for name, kit in SMALLRNA.items() if 'blockers' in kit} == {V4 + ' SE': []}


def mqc_config_command(output):
    return ' '.join(rule_block(output, 'multiqc_config').split())


@pytest.mark.parametrize('blockers, shown', [(None, '--blockers "none" '), (['tRNA/YRNA', 'blood miRNA'], '--blockers "tRNA/YRNA, blood miRNA" ')])
def test_v4_blockers_reach_the_report_header(tmp_path, blockers, shown):
    def update(settings):
        settings['libprepkit'] = V4
        if blockers is not None:
            settings['blockers'] = blockers
    returncode, output = dry_run(tmp_path, 'unitas', update, printshell=True, target='multiqc_report')
    assert returncode == 0, output
    assert shown in mqc_config_command(output)


def test_kits_without_blockers_have_no_header_row(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas', printshell=True, target='multiqc_report')  # fixture: NEXTflex v3
    assert returncode == 0, output
    assert '--blockers' not in mqc_config_command(output)


@pytest.mark.parametrize('kit, blockers, message', [
    (V4, ['YRNA'], 'a list of distinct values from tRNA/YRNA, blood miRNA, custom'),
    (V4, None, 'a list of distinct values'), (V4, 'custom', 'a list of distinct values'),
    (V4, ['custom', 'custom'], 'a list of distinct values'),
    (V3, ['tRNA/YRNA'], 'lists no blockers')])
def test_invalid_blockers_stop(tmp_path, kit, blockers, message):
    returncode, output = dry_run(tmp_path, 'unitas', lambda s: s.update(libprepkit=kit, blockers=blockers))
    assert returncode != 0 and message in output
