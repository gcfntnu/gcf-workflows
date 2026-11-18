"""Small-RNA filtering with cutadapt only (filter.trim.quantifier: cutadapt, the kit default since #350): settings,
commands, the truncated 5' adapter set, the merged report, the workflow route, and read outcomes with cutadapt 5.2 on
synthetic reads.

The read-outcome tests run cutadapt from PATH and skip when it is not version 5.2 (CI installs it).
"""
import copy
import gzip
import importlib.util
import json
from pathlib import Path
import shlex
import shutil
import subprocess

import pytest
import yaml

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / 'smallrna/rules/filter/scripts'
KITS = yaml.safe_load((ROOT / 'libprep.config').read_text())
V3 = 'Bioo Scientific NEXTflex Small RNA-Seq Kit v3'
V4 = 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4'
TRUSEQ = 'TruSeq Small RNA Library Prep Kit'
ADAPTER3 = 'TGGAATTCTCGGGTGCCAAGG'
ADAPTER5 = 'GTTCAGAGTTCTACAGTCCGACGATC'
ADAPTER5_V4 = 'TCTTTCCCTACACGACGCTCTTCCGATCT'  # NEXTflex 5' Adapter v4 (manual, Appendix A), Illumina Read 1 sequence
STP = 'GAATTCCACCACGTTCCCGTGG'  # TruSeq Small RNA Stop Oligo (UniVec NGB00791.1)


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


policy = load('cutadapt_config')
merger = load('merge_cutadapt')


def kit_config(kit, **cutadapt):
    entry = copy.deepcopy(KITS[kit + ' SE'])  # tests change the config; KITS is shared
    config = {'libprepkit': kit, 'adapter': entry['adapter'], 'filter': {'trim': entry['filter']['trim']}}
    if cutadapt:
        config['filter']['trim']['cutadapt'] = dict(config['filter']['trim'].get('cutadapt') or {}, **cutadapt)
    return config


# Settings and commands

def test_kit_presets():
    v3, v4, truseq = (policy.cutadapt_settings(kit_config(kit)) for kit in (V3, V4, TRUSEQ))
    assert policy.filter_args(v3) == '-u 4 -u -4 --action=none --discard-trimmed -m 15 --max-n 5 --max-ee 1'
    assert policy.filter_args(v4) == policy.filter_args(truseq) == '--action=none --discard-trimmed -m 15 --max-n 5 --max-ee 1'
    assert (v3['five_prime_piece'], truseq['five_prime_piece'], v4['five_prime_piece']) == (14, 14, 15)
    assert policy.adapter_args(v3) == policy.adapter_args(v4) == "-a 'TGGAATTCTCGGGTGCCAAGG;min_overlap=5'"
    assert v3['adapter5'] == truseq['adapter5'] == ADAPTER5 and v4['adapter5'] == ADAPTER5_V4
    assert all(KITS[kit + ' SE']['filter']['trim']['quantifier'] == 'cutadapt' for kit in (V3, V4, TRUSEQ))
    assert truseq['reagent_oligos'] == [STP] and v3['reagent_oligos'] == v4['reagent_oligos'] == []
    config = kit_config(TRUSEQ, reagent_oligos=None)
    assert policy.cutadapt_settings(config)['reagent_oligos'] == []                                     # null: none
    del config['filter']['trim']['cutadapt']['reagent_oligos']
    assert policy.cutadapt_settings(config)['reagent_oligos'] == []                                     # unset: none (a custom kit)


def test_settings_are_required():
    with pytest.raises(ValueError, match='none for .*Custom'):
        policy.cutadapt_settings({'libprepkit': 'Custom', 'adapter': ADAPTER3, 'filter': {'trim': {}}})


@pytest.mark.parametrize('change, message', [
    ({'adapter5': None}, 'adapter5'), ({'adapter5': 'GTTCAGAGTTC'}, 'at least 12'), ({'adapter5': 'GTTCAGAGTTCUACAG'}, 'adapter5'),
    ({'random_ends': -1}, 'random_ends'), ({'minimum_length': None}, 'minimum_length'), ({'minimum_length': 0}, 'minimum_length'),
    ({'max_n': 'five'}, 'max_n'), ({'max_ee': 0}, 'max_ee'), ({'five_prime_piece': None}, 'five_prime_piece'),
    ({'five_prime_piece': 13}, 'five_prime_piece'), ({'five_prime_piece': 27}, 'five_prime_piece'), ({'trim_args': '-q 20'}, 'unknown setting'),
    ({'reagent_oligos': STP}, 'reagent_oligos'), ({'reagent_oligos': ['GAATTCCACCACG']}, 'at least 14'),
    ({'reagent_oligos': ['GAAUUCCACCACGUUCC']}, 'reagent_oligos'), ({'reagent_oligos': ''}, 'reagent_oligos'),
    ({'reagent_oligos': False}, 'reagent_oligos'), ({'reagent_oligos': {}}, 'reagent_oligos'),
    ({'reagent_oligos': [STP, STP]}, 'more than once')])
def test_invalid_settings_fail(change, message):
    with pytest.raises(ValueError, match=message):
        policy.cutadapt_settings(kit_config(TRUSEQ, **change))


def test_three_prime_adapter_must_be_a_sequence():
    config = kit_config(TRUSEQ)
    config['adapter'] = 'auto'
    with pytest.raises(ValueError, match="3' adapter"):
        policy.cutadapt_settings(config)


def test_optional_expected_errors():
    assert policy.filter_args(policy.cutadapt_settings(kit_config(TRUSEQ, max_ee=1, max_n=None))) == \
        '--action=none --discard-trimmed -m 15 --max-ee 1'


def test_five_prime_fasta_truncates_the_adapter_and_lists_its_pieces():
    records = policy.five_prime_fasta(ADAPTER5, 14).splitlines()[1::2]
    assert all(r.endswith(';max_error_rate=0.1') for r in records)
    truncated = [r for r in records if ';min_overlap=12;' in r]
    pieces = [r.split(';')[0] for r in records if ';min_overlap=14;' in r]
    assert truncated[0] == 'XGTTCAGAGTTCT;min_overlap=12;max_error_rate=0.1'  # non-internal: only at the read start
    assert truncated[-1] == 'X' + ADAPTER5 + ';min_overlap=12;max_error_rate=0.1' and len(truncated) == len(ADAPTER5) - 12 + 1
    assert 'ACAGTCCGACGATC' in pieces and 'GTTCAGAGTTCTAC' in pieces and all(len(p) >= 14 and p in ADAPTER5 for p in pieces)
    assert len(pieces) == len(set(pieces)) == 91  # every piece of 14-26 nt, prefixes included
    result = subprocess.run(['python', str(SCRIPTS / 'cutadapt_config.py'), 'five-prime-fasta', ADAPTER5, '14'],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0 and result.stdout == policy.five_prime_fasta(ADAPTER5, 14)
    v4_pieces = [r.split(';')[0] for r in policy.five_prime_fasta(ADAPTER5_V4, 15).splitlines()[1::2] if ';min_overlap=15;' in r]
    assert all(len(p) >= 15 and p in ADAPTER5_V4 for p in v4_pieces) and len(v4_pieces) == 15 * 16 // 2


def test_three_prime_artefact_fasta_lists_adapter_primer_and_p7_pieces():
    records = policy.three_prime_artefact_fasta(ADAPTER3).splitlines()
    names, sequences = records[0::2], [r.split(';')[0] for r in records[1::2]]
    assert all(r.endswith(';min_overlap=16;max_error_rate=0.1') for r in records[1::2])
    assert sequences[0] == ADAPTER3 + 'AACTCCAGTCAC' and names[0] == '>three_prime_artefact_adapter_0_33'
    assert 'ATCTCGTATGCCGTCTTCTGCTTG' in sequences and 'CCAAGGAACTCCAGTCAC' in sequences  # P7; adapter end + primer start
    assert all(len(s) >= 16 for s in sequences) and len(sequences) == len(set(sequences)) == 18 * 19 // 2 + 9 * 10 // 2
    result = subprocess.run(['python', str(SCRIPTS / 'cutadapt_config.py'), 'artefact-fasta', ADAPTER3],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0 and result.stdout == policy.three_prime_artefact_fasta(ADAPTER3)


def test_artefact_records_are_named_as_the_merge_counts_them():
    names = [line[1:] for line in policy.artefact_fasta(ADAPTER3, [STP]).splitlines()[0::2]]
    assert (policy.THREE_PRIME_ARTEFACT_RECORD, policy.REAGENT_OLIGO_RECORD) == (merger.THREE_PRIME_ARTEFACT_RECORD, merger.REAGENT_OLIGO_RECORD)
    assert sum(n.startswith(merger.REAGENT_OLIGO_RECORD) for n in names) == 2          # the oligo and its end
    assert all(n.startswith((merger.THREE_PRIME_ARTEFACT_RECORD, merger.REAGENT_OLIGO_RECORD)) for n in names)


def test_artefact_fasta_adds_the_reagent_oligos_as_five_prime_adapters():
    fasta = policy.artefact_fasta(ADAPTER3, [STP])
    assert fasta == (policy.three_prime_artefact_fasta(ADAPTER3) + '>reagent_oligo_1_22\n' + STP + ';min_overlap=14;max_error_rate=0.1\n'
                     '>reagent_oligo_1_end_14\n' + STP[-14:] + ';min_overlap=14;max_error_rate=0.1\n')
    result = subprocess.run(['python', str(SCRIPTS / 'cutadapt_config.py'), 'artefact-fasta', ADAPTER3, STP], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0 and result.stdout == fasta
    result = subprocess.run(['python', str(SCRIPTS / 'cutadapt_config.py'), 'artefact-fasta', ADAPTER3, 'GAATTCC'], capture_output=True, text=True, timeout=30)
    assert result.returncode != 0 and 'at least 14 nt' in result.stderr
    result = subprocess.run(['python', str(SCRIPTS / 'cutadapt_config.py'), 'artefact-fasta', ADAPTER3, STP, STP], capture_output=True, text=True, timeout=30)
    assert result.returncode != 0 and 'more than once' in result.stderr


# The merged report

def report(command, input_reads, output, too_short=0, discard_trimmed=None, with_adapter=0, records=None):
    return {'tag': 'Cutadapt report', 'schema_version': [0, 3], 'cutadapt_version': '5.2', 'command_line_arguments': command,
            'input': {'path1': 'in.fastq', 'path2': None, 'paired': False},
            'read_counts': {'input': input_reads, 'filtered': {'too_short': too_short, 'too_long': None, 'too_many_n': 0,
                                                               'too_many_expected_errors': None, 'casava_filtered': None,
                                                               'discard_trimmed': discard_trimmed, 'discard_untrimmed': None},
                            'output': output, 'reverse_complemented': None, 'read1_with_adapter': with_adapter, 'read2_with_adapter': None},
            'basepair_counts': {'input': input_reads * 50, 'input_read1': input_reads * 50, 'quality_trimmed': None, 'output': output * 20},
            'adapters_read1': records if records is not None else [{'name': command[0]}]}


PATHS = {'adapter': 'a.json', 'filter': 'f.json', 'artefacts': 'x.json'}


def test_merged_report_sums_the_passes():
    adapter, filtered = report(['adapter'], 1000, 1000, with_adapter=990), report(['filter'], 1000, 600, too_short=300, discard_trimmed=100, with_adapter=150)
    records = [{'name': 'three_prime_artefact_adapter_0_33', 'total_matches': 6}, {'name': 'three_prime_artefact_p7_0_24', 'total_matches': 1},
               {'name': 'reagent_oligo_1_22', 'total_matches': 3}]
    artefacts = report(['artefacts'], 600, 590, discard_trimmed=10, with_adapter=10, records=records)
    merged = merger.merge(adapter, filtered, artefacts, 'S1', PATHS)
    counts = merged['read_counts']
    assert counts['input'] == 1000 and counts['output'] == 590 and counts['read1_with_adapter'] == 990
    assert counts['filtered']['too_short'] == 300 and counts['filtered']['discard_trimmed'] == 100  # 5' adapter artefacts
    assert counts['filtered']['gcf_three_prime_artefacts'] == 7 and counts['filtered']['gcf_reagent_oligos'] == 3  # by record
    assert counts['output'] + sum(v or 0 for v in counts['filtered'].values()) == counts['input']
    assert merged['input']['path1'] == 'S1.fastq' and merged['adapters_read1'] == [{'name': 'adapter'}]
    assert merged['tag'].startswith('Cutadapt report, merged') and [p['file'] for p in merged['gcf_passes']] == ['a.json', 'f.json', 'x.json']
    # gcf-bfq's analysis QC also reads these: single-end, and bases out from the last pass
    assert merged['input']['paired'] is False and merged['basepair_counts']['output'] == 590 * 20
    with pytest.raises(ValueError, match='the filter pass read 999 reads, the adapter pass wrote 1000'):
        merger.merge(adapter, report(['filter'], 999, 600), artefacts, 'S1', PATHS)
    with pytest.raises(ValueError, match='the artefacts pass read 600 reads, the filter pass wrote 599'):
        merger.merge(adapter, report(['filter'], 1000, 599), artefacts, 'S1', PATHS)
    with pytest.raises(ValueError, match='other than by discard_trimmed'):
        merger.merge(adapter, filtered, report(['artefacts'], 600, 590, too_short=10), 'S1', PATHS)
    with pytest.raises(ValueError, match='records of unknown kind: other_1'):
        merger.merge(adapter, filtered, report(['artefacts'], 600, 590, discard_trimmed=10, records=records + [{'name': 'other_1', 'total_matches': 0}]),
                     'S1', PATHS)
    with pytest.raises(ValueError, match='discarded 10 reads, but its records matched 9'):
        merger.merge(adapter, filtered, report(['artefacts'], 600, 590, discard_trimmed=10, records=records[1:] + [dict(records[0], total_matches=5)]),
                     'S1', PATHS)


# The workflow route

def use_cutadapt(kit):
    def update(settings):
        settings['libprepkit'] = kit
        settings.setdefault('filter', {}).setdefault('trim', {})['quantifier'] = 'cutadapt'
    return update


def test_route_filters_with_cutadapt_only(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,mirge3', use_cutadapt(V3), printshell=True, target='multiqc_report')
    assert returncode == 0, output
    assert 'rule smallrna_fastp:' not in output and 'rule smallrna_trim_fastp:' not in output
    adapter = rule_block(output, 'smallrna_cutadapt_adapter')
    assert "-a 'TGGAATTCTCGGGTGCCAAGG;min_overlap=5'" in adapter
    filtered = rule_block(output, 'smallrna_cutadapt_filter')
    assert '-g file:data/tmp/smallrna/filter/cutadapt/five_prime_adapters.fa -u 4 -u -4 --action=none --discard-trimmed -m 15 --max-n 5 --max-ee 1' in filtered
    assert 'five-prime-fasta GTTCAGAGTTCTACAGTCCGACGATC 14 ' in rule_block(output, 'smallrna_five_prime_adapters')
    artefacts = rule_block(output, 'smallrna_cutadapt_artefacts')
    assert 'input: data/tmp/smallrna/filter/cutadapt_filtered/SR0' in artefacts
    assert '-g file:data/tmp/smallrna/filter/cutadapt/artefacts.fa --action=none --discard-trimmed' in artefacts
    command = rule_block(output, 'smallrna_artefact_fasta').split('artefact-fasta', 1)[1].split('>', 1)[0].split()
    assert command == ['TGGAATTCTCGGGTGCCAAGG']                                                        # no reagent oligos (NEXTflex v3)
    assert 'bfq/provenance/cutadapt/SR01.artefacts.json' in output
    assert 'input: data/tmp/smallrna/filter/cutadapt/SR0' in rule_block(output, 'unitas').split('output:', 1)[0]  # any sample's job
    assert 'bfq/logs/SR01/SR01.cutadapt.json' in output and 'bfq/provenance/cutadapt/SR01.filter.json' in output
    assert 'SR01.fastp.json' not in output
    multiqc = rule_block(output, 'multiqc_report')
    assert '-m gcf_smallrna_cutadapt' in multiqc and '-m cutadapt ' not in multiqc  # gcf-multiqc's small-RNA module (#358)


def test_route_sets_the_minimum_length_of_the_quantifiers(tmp_path):
    returncode, output = dry_run(tmp_path, 'mirge3', use_cutadapt(TRUSEQ), printshell=True)
    assert returncode == 0, output
    assert '--minimum-length 15 ' in rule_block(output, 'mirge3_quant')
    assert '-u 4' not in rule_block(output, 'smallrna_cutadapt_filter')
    assert 'artefact-fasta TGGAATTCTCGGGTGCCAAGG ' + STP + ' > ' in rule_block(output, 'smallrna_artefact_fasta')  # the Stop Oligo


@pytest.mark.parametrize('kit', [V3, V4, TRUSEQ])
def test_cutadapt_is_the_default(tmp_path, kit):
    returncode, output = dry_run(tmp_path, 'unitas', lambda s: s.update(libprepkit=kit), printshell=True, target='multiqc_report')
    assert returncode == 0, output
    assert 'rule smallrna_cutadapt_filter:' in output and 'rule smallrna_fastp:' not in output
    adapter5 = ADAPTER5_V4 + ' 15 ' if kit == V4 else ADAPTER5 + ' 14 '
    assert 'five-prime-fasta ' + adapter5 in rule_block(output, 'smallrna_five_prime_adapters')


def test_projects_with_fastp_settings_keep_the_fastp_route(tmp_path):
    """Projects initialized before #350 carry the kit's fastp route and settings in config.yaml."""
    def update(settings):
        settings.setdefault('filter', {})['trim'] = {'quantifier': 'fastp', 'fastp': {'args': '-f 4 -t 4 -l 16 -q 15 -u 40 -n 5 -e 0'}}
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    assert 'rule smallrna_fastp:' in output and 'rule smallrna_cutadapt_filter:' not in output
    assert '--minimum-length 16 ' in rule_block(output, 'mirge3_quant')


@pytest.mark.parametrize('quantifier', [None, 'cutadapt', 'skip'])
def test_fastp_settings_with_another_route_stop(tmp_path, quantifier):
    def update(settings):
        trim = settings.setdefault('filter', {}).setdefault('trim', {})
        trim['fastp'] = {'args': '-l 15 -q 15'}
        if quantifier:
            trim['quantifier'] = quantifier
    returncode, output = dry_run(tmp_path, 'unitas', update)
    assert returncode != 0
    assert "filter.trim.fastp.args is set, but filter.trim.quantifier is '{}', so fastp does not run".format(quantifier or 'cutadapt') in output


# Read outcomes with cutadapt 5.2

MIR22_5P, MIR139_5P, MIR21_5P = 'AGTTCTTCAGTGGCAAGCTTTA', 'TCTACAGTGCACGTGTCTCCAG', 'TAGCTTATCAGACTGATGTTGA'
ARTEFACTS = ['GAGTTCTACAGTCCGACGA', 'TACAGTCCGACGATCGCCA', 'AGTTCTACAGTCCGACGATCGGGA', 'GAGTTCTACAGTCTGACGA',
             'GTTCTACAGTCCGACGATCTGGGGTTCTACAGTCCGACGATC', 'AAACAGTCCGACGATCA', 'AATACAGTCCGACGAT']  # last two: after junk bases


def cutadapt_52():
    path = shutil.which('cutadapt')
    if not path:
        return None
    version = subprocess.run([path, '--version'], capture_output=True, text=True).stdout.strip()
    return path if version == '5.2' else None


def run_kit(tmp_path, kit, inserts):
    """Reads [4N]insert[4N]adapter (4N for random_ends 4) through the three passes as the rules run them."""
    cutadapt = cutadapt_52()
    if not cutadapt:
        pytest.skip('cutadapt 5.2 not on PATH')
    settings = policy.cutadapt_settings(kit_config(kit))
    ends = 'ACGT'[:settings['random_ends']]
    raw = tmp_path / 'raw.fastq.gz'
    with gzip.open(raw, 'wt') as fh:
        for i, insert in enumerate(inserts):
            read = (ends + insert + ends[::-1] + ADAPTER3 + 'ATCTCGTATG')[:75]
            fh.write('@r{}\n{}\n+\n{}\n'.format(i, read, 'I' * len(read)))
    (tmp_path / 'five.fa').write_text(policy.five_prime_fasta(settings['adapter5'], settings['five_prime_piece']))
    (tmp_path / 'three.fa').write_text(policy.artefact_fasta(settings['adapter3'], settings['reagent_oligos']))
    first = [cutadapt, *shlex.split(policy.adapter_args(settings)), '--json', str(tmp_path / 'a.json'), '-o', str(tmp_path / 'a.fastq'), str(raw)]
    second = [cutadapt, '-g', 'file:' + str(tmp_path / 'five.fa'), *shlex.split(policy.filter_args(settings)),
              '--json', str(tmp_path / 'f.json'), '-o', str(tmp_path / 'f.fastq'), str(tmp_path / 'a.fastq')]
    third = [cutadapt, '-g', 'file:' + str(tmp_path / 'three.fa'), *shlex.split(policy.artefact_args()),
             '--json', str(tmp_path / 'x.json'), '-o', str(tmp_path / 'x.fastq'), str(tmp_path / 'f.fastq')]
    for command in (first, second, third):
        result = subprocess.run(command, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stderr
    kept = (tmp_path / 'x.fastq').read_text().splitlines()[1::4]
    reports = {name: json.loads((tmp_path / PATHS[name]).read_text()) for name in PATHS}
    return kept, merger.merge(reports['adapter'], reports['filter'], reports['artefacts'], 'S1', PATHS)


@pytest.mark.parametrize('kit', [V3, TRUSEQ])
def test_artefacts_go_and_real_mirnas_stay(tmp_path, kit):
    """miR-22-5p and miR-139-5p begin like part of the 5' adapter (one mismatch over 11 nt): kept with the 12-nt overlap."""
    inserts = [MIR22_5P] * 3 + [MIR139_5P] * 2 + [MIR21_5P] + ARTEFACTS + ['']  # '' is an adapter dimer
    kept, merged = run_kit(tmp_path, kit, inserts)
    assert sorted(kept) == sorted([MIR22_5P] * 3 + [MIR139_5P] * 2 + [MIR21_5P])
    counts = merged['read_counts']
    assert counts['input'] == len(inserts) and counts['output'] == 6
    assert counts['filtered']['discard_trimmed'] == len(ARTEFACTS) and counts['filtered']['too_short'] == 1
    assert counts['filtered']['gcf_three_prime_artefacts'] == 0


# NEXTflex v4: its 5' adapter is the Illumina Read 1 sequence; hsa-miR-6865-3p matches its 14-nt pieces with one error
MIR6865_3P = 'ACACCCTCTTTCCCTACCGCC'
V4_ARTEFACTS = ['CCTACACGACGCTCTTCCG', 'TCTTTCCCTACACGACGC', 'AAAACACGACGCTCTTCCGATC', 'CTACACGACGCTCTTCCGATCTGGA']


def test_v4_artefacts_go_and_mir_6865_3p_stays(tmp_path):
    inserts = [MIR6865_3P] * 2 + [MIR21_5P, MIR22_5P] + V4_ARTEFACTS + ['']
    kept, merged = run_kit(tmp_path, V4, inserts)
    assert sorted(kept) == sorted([MIR6865_3P] * 2 + [MIR21_5P, MIR22_5P])
    assert merged['read_counts']['filtered']['discard_trimmed'] == len(V4_ARTEFACTS)


@pytest.mark.parametrize('piece, kept', [(14, False), (15, True)])
def test_v4_needs_15_nt_pieces(tmp_path, piece, kept):
    cutadapt = cutadapt_52()
    if not cutadapt:
        pytest.skip('cutadapt 5.2 not on PATH')
    (tmp_path / 'five.fa').write_text(policy.five_prime_fasta(ADAPTER5_V4, piece))
    (tmp_path / 'r.fastq').write_text('@r\n{}\n+\n{}\n'.format(MIR6865_3P, 'I' * len(MIR6865_3P)))
    result = subprocess.run([cutadapt, '-g', 'file:' + str(tmp_path / 'five.fa'), '--action=none', '--discard-trimmed',
                             str(tmp_path / 'r.fastq')], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[1::4] == ([MIR6865_3P] if kept else [])


# The 3' artefact pass (#353): reads seen in public plasma data after the 3' adapter and 5' artefact passes
THREE_PRIME_DISCARDED = ['GGCTGGTCCGATGGTAGTGGGTTATCAGAACTTGAATTCTAGGGTGCCAA',   # RNY4 + adapter with 3 errors
                         'TGGAATTCTCGGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTC',   # adapter dimer, adapter with a deletion
                         'TGGAATTCTCGGGTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTC',   # adapter dimer with a duplication
                         'GAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTGAAAAAAA',    # index primer and P7 only
                         'GGGGGGGGGGATCTCGTATGCCGTCTTCTGCTTGGGG']                  # P7 after other bases
THREE_PRIME_KEPT = [MIR21_5P, MIR22_5P, 'GGCTGGTCCGATGGTAGTGGGTTATCAGAACT',  # RNY4 5' fragment
                    'TTTTTAACTCCAGTCATTTTTT',                                 # 12 nt of primer, flanks differ
                    'AAAACCAAGGAACTCCAGAAAA']                                 # 14 nt of adapter end + primer, flanks differ


def test_three_prime_artefacts_are_discarded(tmp_path):
    cutadapt = cutadapt_52()
    if not cutadapt:
        pytest.skip('cutadapt 5.2 not on PATH')
    reads = THREE_PRIME_KEPT + THREE_PRIME_DISCARDED
    (tmp_path / 'three.fa').write_text(policy.three_prime_artefact_fasta(ADAPTER3))
    (tmp_path / 'r.fastq').write_text(''.join('@r{}\n{}\n+\n{}\n'.format(i, r, 'I' * len(r)) for i, r in enumerate(reads)))
    result = subprocess.run([cutadapt, '-g', 'file:' + str(tmp_path / 'three.fa'), *shlex.split(policy.artefact_args()),
                             '--json', str(tmp_path / 'x.json'), str(tmp_path / 'r.fastq')], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[1::4] == THREE_PRIME_KEPT
    assert json.loads((tmp_path / 'x.json').read_text())['read_counts']['filtered']['discard_trimmed'] == len(THREE_PRIME_DISCARDED)


# TruSeq Stop Oligo reads (#380), as in public TruSeq plasma data after the 3' adapter pass: the oligo and its 5'-truncated
# copies (down to 14 nt), with one mismatch (STP[10] is A) or an extra 3' base, and an oligo-miRNA chimera
# the oligo's end counts from 14 nt, at the read start or after other bases (seen in plasma: AA + STP[6:], a copy with
# TTC deleted; one error allowed, so 13 nt after another base too); 13 nt at the read start or 12 nt after others do not
STP_DISCARDED = [STP, STP[4:], STP[6:], STP[-14:] + 'GCATGCA', STP[:10] + 'T' + STP[11:], STP[5:] + 'A', STP + MIR21_5P,
                 'AA' + STP[6:], 'TTTTTTTTT' + STP[-16:], 'TTTTTTTTT' + STP[-13:]]
STP_KEPT = [MIR21_5P, MIR22_5P, MIR139_5P, STP[-13:] + 'GCATGCAT', 'TTTTTTTTTT' + STP[-12:]]


def test_truseq_stop_oligo_reads_are_discarded(tmp_path):
    inserts = STP_KEPT + STP_DISCARDED
    kept, merged = run_kit(tmp_path, TRUSEQ, inserts)
    assert sorted(kept) == sorted(STP_KEPT)
    filtered = merged['read_counts']['filtered']
    assert filtered['gcf_reagent_oligos'] == len(STP_DISCARDED) and filtered['gcf_three_prime_artefacts'] == 0   # counted apart
    kept, _ = run_kit(tmp_path, V3, [STP[4:]])                                                  # no reagent oligo for NEXTflex
    assert kept == [STP[4:]]


# Where and with how many errors the 5' adapter set matches (review of #351): 12-13 nt only at the read start, 14 nt
# anywhere, 10% errors (one below 20 nt, two from 20 nt).
FIVE_PRIME_KEPT = ['CCCCGTTCAGAGTTCTGGCCCC',          # 12-nt prefix inside the read (next two bases differ)
                   'GTTCACAGTTCAACAGTCCTTTTTTTTTT']   # 19-nt prefix at the read start with two mismatches
FIVE_PRIME_DISCARDED = ['CCCGTTCAGAGTTCTACCCCC',     # 14-nt prefix inside the read
                        'CCCCCCCCGTTCAGAGTTCTA',      # 13-nt prefix at the read end: a 14-nt piece with one deletion
                        'CCCCGTTCAGAGTTCTCCCC',       # 14 nt inside the read with one mismatch
                        'GTTCAGAGTTGTGGGGAAAAAA',     # 12-nt prefix at the read start with one mismatch
                        'TTCAGAGTACAACAGTCCGACG']     # 22 nt at the read start with two mismatches


def test_five_prime_matches_by_position_and_errors(tmp_path):
    cutadapt = cutadapt_52()
    if not cutadapt:
        pytest.skip('cutadapt 5.2 not on PATH')
    reads = FIVE_PRIME_KEPT + FIVE_PRIME_DISCARDED
    (tmp_path / 'five.fa').write_text(policy.five_prime_fasta(ADAPTER5, 14))
    (tmp_path / 'r.fastq').write_text(''.join('@r{}\n{}\n+\n{}\n'.format(i, r, 'I' * len(r)) for i, r in enumerate(reads)))
    result = subprocess.run([cutadapt, '-g', 'file:' + str(tmp_path / 'five.fa'), '--action=none', '--discard-trimmed',
                             str(tmp_path / 'r.fastq')], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[1::4] == FIVE_PRIME_KEPT


def test_three_prime_adapter_with_a_deletion_is_trimmed(tmp_path):
    cutadapt = cutadapt_52()
    if not cutadapt:
        pytest.skip('cutadapt 5.2 not on PATH')
    settings = policy.cutadapt_settings(kit_config(TRUSEQ))
    read = MIR21_5P + ADAPTER3.replace('TCGGG', 'TCGG', 1) + 'ATCTCG'  # one base of the adapter missing
    (tmp_path / 'r.fastq').write_text('@r\n{}\n+\n{}\n'.format(read, 'I' * len(read)))
    result = subprocess.run([cutadapt, *shlex.split(policy.adapter_args(settings)), str(tmp_path / 'r.fastq')],
                            capture_output=True, text=True, timeout=60)
    assert result.stdout.splitlines()[1] == MIR21_5P
