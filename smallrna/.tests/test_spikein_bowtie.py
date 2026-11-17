"""Spike-in counts per spike-in with bowtie (filter.spikein.quantifier: bowtie, #336): the reference, count tables and
matrix; the route's dry runs, removal of spike-in reads, project FASTA sets and load errors; and read outcomes with
bowtie 1.x on synthetic reads (skipped when bowtie is not on PATH)."""
import importlib.util
import io
import re
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/filter/scripts/spikein_bowtie.py'
spec = importlib.util.spec_from_file_location('spikein_bowtie', SCRIPT)
spikein = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spikein)

CAL01, CAL02 = 'GTCCCACTCCGTAGATCTGTTC', 'GATGTAACGAGTTGGAATGCAA'  # smallrna_calibrators Cal01, Cal02 (DNA)


# Reference, counts and matrix

def test_reference_combines_sets_as_dna(tmp_path):
    (tmp_path / 'a.fa').write_text('>spikeIn|Cal01 spikeIn 01\nGUCCCACUCCGUAGAUCUGUUC\n')
    (tmp_path / 'b.fa').write_text('>ERCC-00002 x\nacgt\nACGT\n')
    assert spikein.combine([tmp_path / 'a.fa', tmp_path / 'b.fa']) == [('spikeIn|Cal01', CAL01), ('ERCC-00002', 'ACGTACGT')]


@pytest.mark.parametrize('second, message', [('>S1\nACGT\n', 'in both'), (None, 'listed twice')])
def test_reference_names_are_unique(tmp_path, second, message):
    (tmp_path / 'a.fa').write_text('>S1\nACGT\n' + ('>S1 again\nTTTT\n' if second is None else ''))
    paths = [tmp_path / 'a.fa']
    if second:
        (tmp_path / 'b.fa').write_text(second)
        paths.append(tmp_path / 'b.fa')
    with pytest.raises(ValueError, match=message):
        spikein.combine(paths)


def test_counts_list_every_spikein_and_ambiguous_reads():
    hits = io.StringIO('r1\t+\tS1\t0\tACGT\tIIII\t0\n'
                       'r2\t+\tS1\t2\tGT\tII\t0\n'
                       'r3\t+\tS1\t0\tAC\tII\t1\nr3\t+\tS3\t5\tAC\tII\t1\n'   # two spike-ins: ambiguous
                       'r4\t+\tS3\t1\tCC\tII\t1\nr4\t+\tS3\t9\tCC\tII\t1\n')  # twice in one spike-in: counted
    assert spikein.count(['S1', 'S2', 'S3'], hits) == [('S1', 2), ('S2', 0), ('S3', 1), ('__ambiguous', 1)]


@pytest.mark.parametrize('line, message', [('r1\t-\tS1\t0\tACGT\n', 'forward-strand'), ('r1\t+\tS9\t0\tACGT\n', 'not in the spike-in')])
def test_counts_reject_unexpected_alignments(line, message):
    with pytest.raises(ValueError, match=message):
        spikein.count(['S1'], io.StringIO(line))


def test_matrix_puts_samples_side_by_side(tmp_path):
    (tmp_path / 'SR01.tsv').write_text('spikein\tcount\nS1\t3\nS2\t0\n__ambiguous\t1\n')
    (tmp_path / 'SR02.tsv').write_text('spikein\tcount\nS1\t0\nS2\t5\n__ambiguous\t0\n')
    assert spikein.matrix([str(tmp_path / 'SR01.tsv'), str(tmp_path / 'SR02.tsv')]) == \
        'spikein\tSR01\tSR02\nS1\t3\t0\nS2\t0\t5\n__ambiguous\t1\t0\n'
    (tmp_path / 'SR03.tsv').write_text('spikein\tcount\nS2\t0\nS1\t3\n__ambiguous\t1\n')
    with pytest.raises(ValueError, match='SR03'):
        spikein.matrix([str(tmp_path / 'SR01.tsv'), str(tmp_path / 'SR03.tsv')])


# The route in the workflow

def spikein_settings(**spikein_config):
    return lambda settings: settings.setdefault('filter', {}).setdefault('spikein', {}).update(spikein_config)


def test_route_counts_every_spikein(tmp_path):
    """filter.spikein.ref alone turns spike-ins on, with the bowtie route as the default quantifier."""
    returncode, output = dry_run(tmp_path, 'unitas', spikein_settings(ref='ERCC92,smallrna_calibrators'),
                                 printshell=True, target='multiqc_report')
    assert returncode == 0, output
    fasta = rule_block(output, 'spikein_bowtie_fasta')
    assert 'data/ext/spikein/fasta/ERCC92.fa data/ext/spikein/fasta/smallrna_calibrators.fa' in fasta
    bowtie = rule_block(output, 'spikein_bowtie')
    assert 'input: data/tmp/smallrna/filter/cutadapt/SR0' in bowtie
    assert re.search(r'bowtie -q -v 0 --norc -a -p \d+ --un data/tmp/smallrna/filter/spikein/bowtie/SR0', bowtie)
    assert 'bowtie-build data/tmp/smallrna/filter/spikein/spikein.fa data/tmp/smallrna/filter/spikein/index/spikein' in output
    assert 'data/tmp/smallrna/bfq/exprs/spikein_counts.tsv' in rule_block(output, 'bfq_level2_spikein')
    assert 'rule spikein_unitas:' not in output and 'rule unitas_extra_reference:' not in output  # UNITAS: no spike-ins
    mqc = rule_block(output, 'bfq_level2_spikein_mqc')
    assert 'output: data/tmp/smallrna/bfq/figs/spikein_mqc.yaml' in mqc
    assert re.search(r'spikein_mqc.py bowtie --counts data/tmp/smallrna/filter/spikein/bowtie/SR01.tsv .* --logs '
                     r'data/tmp/smallrna/filter/spikein/bowtie/SR01.log', mqc)


def test_bowtie_removes_spikein_reads_for_every_consumer_and_none_estimates_them(tmp_path):
    """bowtie is a prefilter: the quantifiers and miRTrace read the reads without spike-ins, UNITAS gets no spike-in
    reference and miRge3 no spike-in index."""
    update = spikein_settings(quantifier='bowtie', ref='smallrna_calibrators')
    returncode, output = dry_run(tmp_path, 'unitas,mirge3', update, printshell=True, target='multiqc_report')
    assert returncode == 0, output
    assert 'input: data/tmp/smallrna/filter/cutadapt/SR0' in rule_block(output, 'spikein_bowtie')
    for rule in ('unitas', 'qc_mirtrace'):
        assert 'data/tmp/smallrna/filter/spikein/bowtie/SR0' in rule_block(output, rule).split('output:', 1)[0], rule
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)  # quant_all: miRge3 is not delivered
    assert returncode == 0, output
    assert 'ln -sf "$PWD"/data/tmp/smallrna/filter/spikein/bowtie/SR01_R1.fastq' in rule_block(output, 'mirge3_quant')
    assert '--spikeIn' not in rule_block(output, 'mirge3_quant')
    library = [line for line in rule_block(output, 'mirge3_project_library').splitlines() if 'mirge3_project_library.py' in line][0]
    assert library.rstrip().endswith('--record data/tmp/smallrna/quant/mirge3/spikein.tsv')  # no spike-in FASTA


def test_native_spikeins_reach_only_quantifiers_that_estimate_them(tmp_path):
    """native: UNITAS and miRge3 read the reads with the spike-ins and estimate them; mirflowz, miRTrace and STAR read
    the bowtie prefilter's reads without spike-ins; bowtie's counts per spike-in are delivered too."""
    update = spikein_settings(quantifier='native', ref='smallrna_calibrators')
    trimmed, without = 'data/tmp/smallrna/filter/cutadapt/SR0', 'data/tmp/smallrna/filter/spikein/bowtie/SR0'
    returncode, output = dry_run(tmp_path, 'unitas', update, printshell=True, target='multiqc_report')
    assert returncode == 0, output
    assert trimmed in rule_block(output, 'unitas').split('output:', 1)[0] and '-refseq' in rule_block(output, 'unitas')
    for rule in ('qc_mirtrace', 'star_align'):
        assert without in rule_block(output, rule).split('output:', 1)[0], rule
    assert 'data/tmp/smallrna/bfq/exprs/spikein_counts.tsv' in rule_block(output, 'bfq_level2_spikein')
    assert 'spikein_mqc.py unitas ' in rule_block(output, 'bfq_level2_spikein_mqc')  # General Statistics from UNITAS
    returncode, output = dry_run(tmp_path, 'mirge3,mirflowz', update, printshell=True)
    assert returncode == 0, output
    assert 'ln -sf "$PWD"/data/tmp/smallrna/filter/cutadapt/SR01_R1.fastq' in rule_block(output, 'mirge3_quant')
    assert '--spikeIn' in rule_block(output, 'mirge3_quant')
    assert without in rule_block(output, 'mirflowz_quant') and trimmed not in rule_block(output, 'mirflowz_quant')


def test_project_fasta_reaches_every_route_in_dna_letters(tmp_path):
    """A project FASTA (here RNA letters, a space in its name) is used as its DNA copy by the bowtie route, UNITAS and
    miRge3; the copying command quotes the path and writes DNA."""
    (tmp_path / 'my spikeins.fa').write_text('>Syn1 synthetic\nACGUACGUACGUACGUACGU\n')
    copy = 'data/tmp/smallrna/filter/spikein/sets/custom1.fa'
    for quantifier, method, rule, target in (('bowtie', 'unitas', 'spikein_bowtie_fasta', 'multiqc_report'),
                                             ('native', 'unitas', 'spikein_unitas', 'multiqc_report'),
                                             ('native', 'mirge3', 'mirge3_project_library', 'quant_all')):
        returncode, output = dry_run(tmp_path, method, spikein_settings(quantifier=quantifier, ref='smallrna_calibrators,my spikeins.fa'),
                                     printshell=True, target=target)
        assert returncode == 0, output
        block = rule_block(output, rule)
        assert 'data/ext/spikein/fasta/smallrna_calibrators.fa' in block and copy in block and 'my spikeins' not in block, rule
    command = [line for line in rule_block(output, 'spikein_project_fasta').splitlines() if 'spikein_bowtie.py' in line][0]
    command = command.split('Shell command:', 1)[-1].strip()
    assert "dna 'my spikeins.fa' > " in command
    result = subprocess.run(command.replace(' > ' + copy, ''), shell=True, cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0 and result.stdout == '>Syn1 synthetic\nACGTACGTACGTACGTACGT\n', result.stderr


def test_spikein_mqc_table_for_each_route(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas', spikein_settings(quantifier='native', ref='ERCC92'), printshell=True,
                                 target='multiqc_report')
    assert returncode == 0, output
    assert re.search(r'spikein_mqc.py unitas data/tmp/smallrna/quant/unitas/unitas/SR01/unitas.annotation_summary.txt',
                     rule_block(output, 'bfq_level2_spikein_mqc'))
    returncode, output = dry_run(tmp_path / 'none', 'unitas', spikein_settings(ref=None), target='multiqc_report')  # fresh config
    assert returncode == 0, output
    assert 'spikein_mqc' not in output and 'spikein_bowtie' not in output


@pytest.mark.parametrize('spikein_config, message', [
    ({'quantifier': 'unitas', 'ref': 'ERCC92'}, 'filter.spikein.quantifier: unitas is now called native'),
    ({'quantifier': 'bowtie', 'ref': 'ERCC92', 'filter': True}, 'filter.spikein.filter is no longer used'),
    ({'quantifier': 'native', 'ref': 'ERCC92', 'filter': False}, 'filter.spikein.filter is no longer used'),
    ({'quantifier': 'skip', 'ref': 'ERCC92'}, 'filter.spikein.quantifier: skip is no longer used'),
    ({'quantifier': 'skip', 'ref': None}, 'filter.spikein.quantifier: skip is no longer used'),
    ({'quantifier': 'native', 'ref': None}, 'filter.spikein.quantifier is native but filter.spikein.ref lists no spike-in sets'),
    ({'quantifier': 'bowtie2', 'ref': 'ERCC92'}, 'filter.spikein.quantifier: smallrna supports bowtie, native'),
    ({'quantifier': 'bowtie', 'ref': 'ERCC92,ERCC92'}, 'listed twice'),
    ({'quantifier': 'bowtie', 'ref': 'missing.fa'}, "'missing.fa' is neither a gcfdb spike-in set"),
    ({'quantifier': 'bowtie', 'ref': 'bad.fa'}, "record 'Syn1' has bases other than A/C/G/T/U"),
    ({'quantifier': 'bowtie', 'ref': 'empty.fa'}, "record 'Syn2' has no sequence")])
def test_invalid_spikein_settings_stop_at_load(tmp_path, spikein_config, message):
    (tmp_path / 'bad.fa').write_text('>Syn1\nACGTNNACGT\n')
    (tmp_path / 'empty.fa').write_text('>Syn1\nACGT\n>Syn2\n')
    returncode, output = dry_run(tmp_path, 'unitas', spikein_settings(**spikein_config))
    assert returncode != 0 and message in output, output[-2000:]


# MultiQC General Statistics

MQC = importlib.util.spec_from_file_location('spikein_mqc', ROOT / 'smallrna/rules/bfq/scripts/spikein_mqc.py')
spikein_mqc = importlib.util.module_from_spec(MQC)
MQC.loader.exec_module(spikein_mqc)


def test_mqc_totals_from_bowtie_tables_and_logs(tmp_path):
    (tmp_path / 'SR01.tsv').write_text('spikein\tcount\nS1\t20\nS2\t3\n__ambiguous\t1\n')
    (tmp_path / 'SR01.log').write_text('# reads processed: 200\n# reads with at least one reported alignment: 24 (12.00%)\n')
    totals = spikein_mqc.bowtie_totals([str(tmp_path / 'SR01.tsv')], [str(tmp_path / 'SR01.log')])
    assert totals == {'SR01': (24, 200)}
    text = spikein_mqc.mqc('bowtie', totals)
    assert 'plot_type: generalstats' in text and '"SR01": {"spikein_reads": 24, "spikein_percent": 12.0}' in text
    import yaml
    parsed = yaml.safe_load(text)
    assert parsed['data'] == {'SR01': {'spikein_reads': 24, 'spikein_percent': 12.0}}
    assert [list(c)[0] for c in parsed['pconfig']] == ['spikein_reads', 'spikein_percent']


def test_mqc_totals_from_unitas_summaries(tmp_path):
    summary = tmp_path / 'SR02' / 'unitas.annotation_summary.txt'
    summary.parent.mkdir()
    summary.write_text('low_complexity\t0\nmiRNA\t150\n   miRNA:mus_musculus\t150\nSpikein\t27\n'
                       '   spikeIn::Cal01 spikeIn 01\t12\nno annotation\t20.5\n')
    assert spikein_mqc.unitas_totals([str(summary)]) == {'SR02': (27.0, 197.5)}


# Read outcomes with bowtie 1.x

def bowtie():
    path, build = shutil.which('bowtie'), shutil.which('bowtie-build')
    return (path, build) if path and build else None


def test_bowtie_counts_exact_forward_pieces_only(tmp_path):
    tools = bowtie()
    if not tools:
        pytest.skip('bowtie and bowtie-build not on PATH')
    (tmp_path / 'cal.fa').write_text('>spikeIn|Cal01 spikeIn 01\n{}\n>spikeIn|Cal02 spikeIn 2\n{}\n'.format(
        CAL01.replace('T', 'U'), CAL02.replace('T', 'U')))
    (tmp_path / 'own.fa').write_text('>Syn1\nAAAACCCCGGGGTTTTACGTAC\n>Syn2\nTTTTCCCCGGGGTTTTACGTAC\n'  # share 18 nt
                                     '>Rep1\nGATCGATCAAGGCCTTGGGGGGATCGATCAAGGCCTTGG\n'          # motif twice
                                     '>Rep2\nTTTTTGATCGATCAAGGCCTTGGCCCCC\n')                    # motif once
    reference = tmp_path / 'spikein.fa'
    reference.write_text(''.join('>{}\n{}\n'.format(n, s) for n, s in spikein.combine([tmp_path / 'cal.fa', tmp_path / 'own.fa'])))
    rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    reads = {'full': CAL01, 'piece': CAL01[2:20], 'cal02': CAL02, 'mismatch': CAL02[:10] + 'C' + CAL02[11:],
             'reverse': rc(CAL01), 'shared': 'CCCCGGGGTTTTACGTAC', 'mirna': 'TAGCTTATCAGACTGATGTTGA',
             'longer': CAL01 + 'AAAA', 'repeated': 'GATCGATCAAGGCCTTGG'}  # twice in Rep1, once in Rep2: ambiguous
    fastq = tmp_path / 'reads.fastq'
    fastq.write_text(''.join('@{}\n{}\n+\n{}\n'.format(name, seq, 'I' * len(seq)) for name, seq in reads.items()))
    subprocess.run([tools[1], str(reference), str(tmp_path / 'idx')], check=True, capture_output=True, timeout=120)
    hits = tmp_path / 'reads.hits'
    with open(hits, 'w') as out:
        result = subprocess.run([tools[0], '-q', '-v', '0', '--norc', '-a', '-p', '1', '--un', str(tmp_path / 'un.fastq'),
                                 str(tmp_path / 'idx'), str(fastq)], stdout=out, stderr=subprocess.PIPE, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    table = subprocess.run([sys.executable, str(SCRIPT), 'counts', '--fasta', str(reference), '--hits', str(hits)],
                           capture_output=True, text=True, timeout=60)
    assert table.returncode == 0, table.stderr
    assert table.stdout == ('spikein\tcount\nspikeIn|Cal01\t2\nspikeIn|Cal02\t1\nSyn1\t0\nSyn2\t0\nRep1\t0\nRep2\t0\n'
                            '__ambiguous\t2\n')
    unmatched = (tmp_path / 'un.fastq').read_text().splitlines()[0::4]
    assert sorted(unmatched) == ['@longer', '@mirna', '@mismatch', '@reverse']
