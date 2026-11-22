"""Tissue mixtures (#451): the design, the prediction of a mixture from its subject's pure samples, the fitted content
ratio, the expected changes between fractions, and prepare and evaluate of a mixture dataset, on hand-made counts."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / 'smallrna/benchmark'
sys.path.insert(0, str(BENCH))
import benchmark  # noqa: E402
import de  # noqa: E402
import mixture  # noqa: E402
from test_benchmark import unitas_dump, write, write_cutadapt_reports, write_records  # noqa: E402

FRACTIONS = (0, 0.25, 0.5, 0.75, 1)
# miRNA: (proportion in the pure first tissue, in the pure second), before scaling to 1
PROFILE = {'liver_only': (0.60, 0.0), 'spleen_only': (0.0, 0.30), 'both': (0.20, 0.20), 'more_spleen': (0.05, 0.40),
           'more_liver': (0.149, 0.099), 'rare': (0.001, 0.001)}
READS = 200000


def samples(subjects=('A', 'B', 'C'), replicate=True, pool=True):
    rows = [{'sample': '{}_{:03d}'.format(s, round(100 * f)), 'subject': s, 'fraction': f} for s in subjects for f in FRACTIONS]
    if replicate:
        rows.append({'sample': 'A_050_r2', 'subject': 'A', 'fraction': 0.5, 'replicate': 2})
    if pool:
        rows.append({'sample': 'pool_000', 'subject': 'pool', 'fraction': 0})
    return rows


def model_counts(k, design, change=None):
    """{sample: {miRNA: reads}} that follow the model exactly (rounded): a mixture's proportions are w(fraction) of the
    second tissue's and 1 - w of the first's. `change(sample, fraction, counts)` edits a sample's counts."""
    first = {m: a / sum(a for a, _ in PROFILE.values()) for m, (a, _) in PROFILE.items()}
    second = {m: b / sum(b for _, b in PROFILE.values()) for m, (_, b) in PROFILE.items()}
    counts = {}
    for row in design:
        w = mixture.weight(row['fraction'], k[row['subject']] if isinstance(k, dict) else k)
        counts[row['sample']] = {m: float(round(READS * (w * second[m] + (1 - w) * first[m]))) for m in PROFILE}
        if change:
            change(row['sample'], row['fraction'], counts[row['sample']])
    return counts


def scored(k=0.5, change=None, rows=None):
    rows = rows or samples()
    design = mixture.Design(rows)
    counts = mixture.Counts(model_counts(k, rows, change), design.samples, sorted(PROFILE))
    return design, counts, mixture.Scores(counts, design)


def test_design_gives_each_subject_its_series_and_keeps_replicates_and_pools_apart():
    design = mixture.Design(samples())
    assert sorted(design.series) == ['A', 'B', 'C'] and design.fractions == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert design.series['A'] == [(0.0, 'A_000'), (0.25, 'A_025'), (0.5, 'A_050'), (0.75, 'A_075'), (1.0, 'A_100')]  # the first replicate
    assert design.replicates == [('A_050', 'A_050_r2')] and design.unused == ['pool_000']  # a pool has no pure second tissue
    assert design.at(0.5) == [('A', 'A_050'), ('B', 'B_050'), ('C', 'C_050')]
    assert mixture.contrasts(design) == [(1.0, 0.0, False), (0.5, 0.25, True), (0.75, 0.25, True), (0.75, 0.5, True)]
    assert [mixture.label(f) for f in (0, 0.25, 1)] == ['0', '25', '100']
    for change, message in ((lambda rows: rows[0].update(fraction=1.5), 'fraction must be a number from 0 to 1'),
                            (lambda rows: rows[0].update(fraction='half'), 'fraction must be a number from 0 to 1'),
                            (lambda rows: rows[0].pop('subject'), 'needs a subject'),
                            (lambda rows: rows[-2].update(replicate=1), 'two samples at fraction 0.5 with the same replicate number')):
        rows = samples()
        change(rows)
        with pytest.raises(ValueError, match=message):
            mixture.Design(rows)
    with pytest.raises(ValueError, match='at least two subjects'):
        mixture.Design(samples(subjects=('A',)))


def test_counts_that_follow_the_model_are_predicted_and_give_back_the_content_ratio():
    design, counts, scores = scored(k=0.5)
    assert scores.k == pytest.approx(0.5, rel=0.01) and scores.spread < 0.01
    assert mixture.weight(0.5, 0.5) == pytest.approx(1 / 3)  # half the RNA from a tissue with half the miRNA: a third of the reads
    values = scores.values()
    assert values['linearity_error'] < 0.01 and values['within_1_5fold'] == 1.0 and values['monotonic'] == 1.0
    assert values['replicate_error'] == 0.0 and values['content_ratio'] == pytest.approx(0.5, abs=0.01)
    assert int((~mixture.np.isnan(scores.residual)).sum()) == 6 * 9  # every miRNA in each of the 3 x 3 mixtures has 10 predicted reads
    shallow = mixture.Scores(mixture.Counts({s: {m: c // 100 for m, c in row.items()} for s, row in model_counts(0.5, samples()).items()},
                                            design.samples, sorted(PROFILE)), design)
    assert int((~mixture.np.isnan(shallow.residual)).sum()) == 5 * 9  # at a hundredth of the depth `rare` has 2 reads: not scored
    assert [row[0] for row in scores.by_abundance()] == ['10-99', '100-999', '>= 1000']
    # one subject, and a resample of miRNAs (row indices, with repeats), score without error
    assert scores.values(subjects=['B'])['linearity_error'] < 0.01 and scores.values([0, 0, 1], ['A', 'A'])['within_1_5fold'] == 1.0
    assert all(math.isnan(v) for k, v in scores.values([], ['A']).items() if k in mixture.SCORE_KEYS)


def test_each_subject_has_its_own_content_ratio():
    """The tissues of one subject can hold another amount of miRNA per unit of RNA than those of the next: each subject's
    mixtures are predicted with its own ratio, so that difference is not counted as error."""
    rows = samples()
    design = mixture.Design(rows)
    truth = {'A': 0.5, 'B': 2.0, 'C': 1.0, 'pool': 1.0}
    scores = mixture.Scores(mixture.Counts(model_counts(truth, rows), design.samples, sorted(PROFILE)), design)
    assert {s: round(k, 2) for s, k in scores.ks.items()} == {'A': 0.5, 'B': 2.0, 'C': 1.0}
    values = scores.values()
    assert values['linearity_error'] < 0.01 and values['content_ratio'] == pytest.approx(1.0, abs=0.01)
    assert values['content_ratio_spread'] == pytest.approx(1.0, abs=0.02)  # the subjects are a factor 2 from the median


def test_a_constant_bias_changes_the_content_ratio_and_not_the_error():
    """Counting three times too much of the liver's main miRNA in every sample is a per-miRNA bias: the mixtures stay
    weighted sums of the pure samples, with another content ratio."""
    _, _, plain = scored(k=1.0)
    _, _, biased = scored(k=1.0, change=lambda sample, fraction, counts: counts.update(liver_only=3 * counts['liver_only']))
    assert biased.values()['linearity_error'] < 0.02 and biased.values()['monotonic'] == 1.0
    assert biased.k < 0.6 * plain.k  # the liver's reads more than doubled


def test_counting_that_differs_between_samples_shows_as_error():
    """A method that loses half of one miRNA's reads in the mixtures only (not in the pure samples) is not linear."""
    def lose(sample, fraction, counts):
        if 0 < fraction < 1:
            counts['more_spleen'] = counts['more_spleen'] // 2
    _, _, scores = scored(k=1.0, change=lose)
    _, _, plain = scored(k=1.0)
    assert scores.values()['linearity_error'] > 5 * max(plain.values()['linearity_error'], 0.01)
    assert scores.values()['within_1_5fold'] < 0.9

    def swap(sample, fraction, counts):  # a tissue-specific miRNA that rises and falls: not monotonic
        if fraction == 0.5:
            counts['spleen_only'] = counts['spleen_only'] // 10
    assert scored(k=1.0, change=swap)[2].values()['monotonic'] < 1.0


def test_technical_replicates_give_the_noise_floor():
    def noisy(sample, fraction, counts):
        if sample == 'A_050_r2':
            counts['both'] = counts['both'] * 2
    _, _, scores = scored(change=noisy)
    assert 0 < scores.values()['replicate_error'] < 1 and scores.values()['linearity_error'] < 0.01  # the replicate is not in the series
    assert math.isnan(scored(rows=samples(replicate=False))[2].values()['replicate_error'])


def test_expected_changes_between_fractions_come_from_the_pure_samples():
    design, counts, scores = scored(k=0.5)
    expected = mixture.expected_changes(counts, design, scores.ks, 0.75, 0.25)
    w = lambda f: mixture.weight(f, scores.k)
    assert expected['liver_only'] == pytest.approx(math.log2((1 - w(0.75)) / (1 - w(0.25))), abs=0.01)  # falls with the spleen's share
    assert expected['spleen_only'] == pytest.approx(math.log2(w(0.75) / w(0.25)), abs=0.01)
    assert expected['both'] == pytest.approx(0, abs=0.01) and expected['more_spleen'] > 1 > expected['more_liver']
    assert mixture.expected_changes(counts, design, {}, 0.75, 0.25) == {}
    empty = mixture.Counts({}, design.samples, sorted(PROFILE))  # a method without counts: no ratio, no expectation, no error
    assert mixture.Scores(empty, design).ks == {} and mixture.expected_changes(empty, design, dict.fromkeys(design.series, 1.0), 0.75, 0.25) == {}


MATURE = ''.join('>mmu-miR-{0}-5p MIMAT00000{0} Mus musculus miR-{0}-5p\n{1}\n'.format(10 + i, seq) for i, seq in enumerate(
    ['UGGAGUGUGACAAUGGUGUUUG', 'UAGCAGCACGUAAAUAUUGGCG', 'UAGCUUAUCAGACUGAUGUUGA', 'UGAGGUAGUAGGUUGUAUAGUU', 'UUCAAGUAAUCCAGGAUAGGCU',
     'CAUUGCACUUGUCUCGGUCUGA']))
NAMES = dict(zip(sorted(PROFILE), ['miR-{}-5p'.format(10 + i) for i in range(6)]))


def declaration(tmp_path, change=None):
    rows = samples(subjects=('A', 'B'), replicate=True, pool=True)
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/mix'
    write(data / 'mirbase/23.1/mature.fa', MATURE)
    for i, row in enumerate(rows):
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/RUN{}.fastq.gz'.format(i), 'wt') as fh:
            fh.write('@r1\nTGGAGTGTGACAATGGTGTTTGAGATCGGAAG\n+\n{}\n'.format('I' * 32))
        row['runs'] = [{'run': 'RUN{}'.format(i), 'md5': 'unused', 'url': 'unused'}]
    spec = {'name': 'mix', 'kind': 'mixture', 'organism': 'mus_musculus', 'libprepkit': 'NEBNext Multiplex Small RNA Library Prep Kit',
            'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/mix', 'mixture': {'components': ['liver', 'spleen']},
            'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'}, 'de': {}, 'samples': rows}
    if change:
        change(spec)
    return ext, write(tmp_path / 'mix.yaml', yaml.safe_dump(spec, sort_keys=False)), rows


def test_a_mixture_declaration_is_checked_and_gives_its_de_sets(tmp_path):
    _, path, _ = declaration(tmp_path)
    spec = benchmark.load_dataset(path)
    assert benchmark.de_kind(spec) == 'mixture'
    sets = benchmark.mixture_de_sets(spec)
    assert [(s.suffix, s.folder, between) for s, _, _, between in sets] == [
        ('', '', False), ('.50_vs_25', 'contrast/50_vs_25', True), ('.75_vs_25', 'contrast/75_vs_25', True), ('.75_vs_50', 'contrast/75_vs_50', True)]
    pure = sets[0][0]
    assert pure.design == [('A_100', '100', {'subject': 'A'}), ('B_100', '100', {'subject': 'B'}), ('A_000', '0', {'subject': 'A'}),
                           ('B_000', '0', {'subject': 'B'})]  # paired by subject; no replicate and no pool
    assert pure.model == de.model_text('100', '0', (), 10) and pure.samples == '2 100 and 2 0 samples of 2 subjects'
    for change, message in ((lambda s: s['mixture'].update(components=['liver']), 'mixture.components must name the two tissues'),
                            (lambda s: s.pop('mixture'), "a mixture dataset needs 'mixture'"),
                            (lambda s: s['samples'][0].update(fraction=2), 'fraction must be a number from 0 to 1'),
                            (lambda s: s.update(de={'covariates': ['subject']}), '`de` takes no options for this dataset')):
        with pytest.raises(SystemExit, match=message):
            benchmark.load_dataset(declaration(tmp_path, change)[1])
    shipped = benchmark.load_dataset(BENCH / 'datasets/srp697344_liver_spleen_nebnext.yaml')
    design = mixture.Design(shipped['samples'])
    assert sorted(design.series) == ['M1', 'M2', 'M3', 'M4', 'M5'] and design.fractions == [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    assert len(design.replicates) == 5 and design.unused == ['pool_S000'] and len(benchmark.mixture_de_sets(shipped)) == 7


def test_mixture_de_scores_use_each_methods_own_expected_changes():
    features = ['a', 'b', 'c', 'd']
    expected = {'m1': {'.x': {'a': 2.0, 'b': -1.5, 'c': 0.5, 'd': 0.05}}, 'm2': {'.x': {'a': 1.2, 'b': 0.4, 'c': 0.5, 'd': 0.05}}}
    fits = {'m1': {'a': (1.8, 0.001, 0.01), 'b': (1.0, 0.001, 0.01), 'c': (0.4, 0.2, 0.5), 'd': (0.3, 0.001, 0.01)},
            'm2': {'a': (1.0, 0.2, 0.3), 'b': (None, None, None), 'c': (-0.6, 0.001, 0.01), 'd': (0.0, 0.5, 0.9)}}
    results = {'.x': {(m, test, f): fits[m][f] for m in fits for test in de.TESTS for f in features}}
    values, per_contrast = benchmark.mixture_de_scores(results, expected, ['m1', 'm2'])
    assert values['m1']['edger_called'] == 0.5            # a called as expected; b called, but up where a fall is expected
    assert values['m1']['voom_wrong_direction'] == 0.5    # of its calls a and b (d is under the least expected change: not scored)
    assert values['m1']['edger_glm_fc_bias'] == -0.2      # median of -0.2 (a), -2.5 (b) and -0.1 (c)
    assert values['m2']['edger_called'] == 0.0 and values['m2']['edger_wrong_direction'] == 1.0  # its one call, c, against the expectation
    assert values['m1']['edger_called_shared'] == 1.0 and values['m2']['edger_called_shared'] == 0.0  # only a is 2-fold for both
    assert values['m1']['voom_unchanged_called'] == 1.0 and values['m2']['voom_unchanged_called'] == 0.0  # d: hardly any expected change
    assert per_contrast == {'.x': {'m1': (2, 0.5), 'm2': (1, 0.0)}}  # b is untested by m2


def test_prepare_and_evaluate_a_mixture(tmp_path):
    ext, path, rows = declaration(tmp_path)
    project = tmp_path / 'project'
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]
    result = subprocess.run(cmd + ['prepare', '--dataset', str(path), '--project', str(project), '--ext-dir', str(ext), '--lanes'],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and not project.exists()  # --lanes is for biological datasets
    result = subprocess.run(cmd + ['prepare', '--dataset', str(path), '--project', str(project), '--ext-dir', str(ext)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load((project / 'config.yaml').read_text())
    assert config['libprepkit'].startswith('NEBNext') and len(config['samples']) == 12
    names = list(config['samples'])
    quant = project / 'data/tmp/smallrna/quant'
    counts = {'unitas': model_counts(0.5, rows), 'mirge3': model_counts(0.5, rows, lambda s, f, c: c.update(more_spleen=c['more_spleen'] // (2 if 0 < f < 1 else 1)))}
    write(quant / 'unitas/mir_counts.tsv', '\t' + '\t'.join(names) + '\n' + ''.join(
        '{}\t{}\n'.format(NAMES[m], '\t'.join('{:g}'.format(counts['unitas'][s][m]) for s in names)) for m in sorted(PROFILE)))
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,' + ','.join(names) + '\n' + ''.join(
        'mmu-{},{}\n'.format(NAMES[m], ','.join('{:g}'.format(counts['mirge3'][s][m]) for s in names)) for m in sorted(PROFILE)))
    write_records(quant, ['unitas', 'mirge3'], dump=unitas_dump(ext, 'mus_musculus', 'Mouse', {}))
    evaluate = lambda *args: subprocess.run(cmd + ['evaluate', '--project', str(project), '--bootstrap', '50', *args],
                                            capture_output=True, text=True, timeout=120)
    result = evaluate()
    assert result.returncode == 0, result.stderr
    bench = project / 'benchmark'
    table = {line.split('\t')[0]: line.split('\t')[1:] for line in (bench / 'metrics.tsv').read_text().splitlines()}
    assert table['metric'] == ['unitas', 'mirge3'] and float(table['content_ratio'][0]) == pytest.approx(0.5, abs=0.01)
    assert float(table['linearity_error'][0]) < 0.01 < float(table['linearity_error'][1])  # miRge3 here loses reads in the mixtures only
    assert table['within_1_5fold'][0] == '1.0' and table['monotonic'][0] == '1.0' and table['replicate_error'] == ['0.0', '0.0']
    assert 'edger_called' not in table
    report = (bench / 'report.md').read_text()
    assert 'total RNA of liver and spleen of 2 subjects, mixed per subject at 0%, 25%, 50%, 75%, 100% spleen' in report
    assert '1 technical replicate pair(s); outside the series (e.g. a pool of several subjects): pool_000.' in report
    assert '## Linearity' in report and '| A | ' in report and '| >= 1000 |' in report and '## Uncertainty' in report
    assert 'Per subject, `linearity_error` / content ratio' in report and ' / 0.5' in report  # each subject's own ratio
    assert any(line.startswith('linearity_error\tunitas\t') for line in (bench / 'metrics_ci.tsv').read_text().splitlines())
    assert (bench / 'de/levels.tsv').read_text() == 'column\tlevel\ncontrast\t50_vs_25\ncontrast\t75_vs_25\ncontrast\t75_vs_50\n'
    assert (bench / 'de/contrast/75_vs_25/input/design.tsv').read_text() == (
        'sample\tgroup\tsubject\nA_075\t75\tA\nB_075\t75\tB\nA_025\t25\tA\nB_025\t25\tB\n')
    assert (bench / 'de/contrast/75_vs_25/input/model.tsv').read_text() == de.model_text('75', '25', (), 10)
    assert (bench / 'de/input/counts.tsv').read_text().splitlines()[0] == 'row\tfeature\tA_100\tB_100\tA_000\tB_000'
    sizes = (bench / 'de/input/library_sizes.tsv').read_text().splitlines()
    assert sizes[1].split('\t')[1:] == [repr(sum(counts['unitas'][s].values())) for s in ('A_100', 'B_100', 'A_000', 'B_000')]  # placed reads
    result = evaluate('--de')
    assert result.returncode != 0 and 'no results.tsv or record.tsv' in result.stderr

    accession = {m: 'MIMAT00000{}'.format(10 + i) for i, m in enumerate(sorted(PROFILE))}
    w = lambda f: mixture.weight(f, 0.5)

    def outputs(folder, high, low):
        """de.smk's outputs for a set: every method finds the changes the model gives; `rare` is not tested."""
        base = bench / 'de' / folder
        texts = {name: (base / 'input' / name).read_text() for name in de.INPUTS}
        change = {'liver_only': math.log2((1 - w(high)) / (1 - w(low))) if high < 1 else -8.0, 'spleen_only': math.log2(w(high) / w(low)) if low > 0 else 8.0,
                  'both': 0.0, 'more_spleen': 1.2, 'more_liver': -0.3}
        lines = ['\t'.join(de.RESULT_HEADER)]
        for row in ('unitas', 'mirge3'):
            for test in de.TESTS:
                for m in sorted(PROFILE):
                    fit = ('NA', 'NA', 'NA') if m == 'rare' else (change[m], 0.001, 0.01) if abs(change[m]) > 0.25 else (change[m], 0.5, 0.9)
                    lines.append('\t'.join([row, test, accession[m]] + [str(v) for v in fit]))
        text = '\n'.join(lines) + '\n'
        write(base / 'results.tsv', text)
        record = {'md5_' + Path(name).stem: de.md5(t) for name, t in texts.items()}
        record.update(md5_script=de.md5(de.SCRIPT.read_bytes()), md5_results=de.md5(text), R='R version 4.5.3', edgeR='4.10.7', limma='3.68.5',
                      formula='~subject + group', min_mean='10', edger='exact test', edger_glm='GLM test', voom='voom')
        write(base / 'record.tsv', 'field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))

    outputs('', 1.0, 0.0)
    for high, low in ((0.5, 0.25), (0.75, 0.25), (0.75, 0.5)):
        outputs('contrast/{:g}_vs_{:g}'.format(100 * high, 100 * low), high, low)
    result = evaluate('--de')
    assert result.returncode == 0, result.stderr
    table = {line.split('\t')[0]: line.split('\t')[1:] for line in (bench / 'metrics.tsv').read_text().splitlines()}
    assert table['edger_glm_called'][0] == '1.0' and table['edger_glm_wrong_direction'][0] == '0.0'  # UNITAS: every expected change found
    assert table['voom_called_shared'][0] == '1.0' and float(table['edger_fc_bias'][0]) == pytest.approx(0, abs=0.2)
    report = (bench / 'report.md').read_text()
    assert 'pure spleen (100) against pure liver (0) over the 6 miRNAs any method reports' in report
    assert "Library sizes are a method's reads on miRNAs with an accession" in report and 'the exact test takes no covariates' in report
    assert 'edger_unchanged_called' in table
    assert '### Between the mixtures' in report and '| 75 against 25 | 2 75 and 2 25 samples of 2 subjects |' in report
    sha = hashlib.sha256((bench / 'de/contrast/75_vs_25/results.tsv').read_bytes()).hexdigest()
    assert '75_vs_25: design `~subject + group` (edgeR 4.10.7, limma 3.68.5); `de/contrast/75_vs_25/results.tsv` SHA-256 `{}`.'.format(sha) in report
    record = json.loads((bench / 'evaluation.json').read_text())
    assert record['bootstrap']['resamples'] == 50
    out = tmp_path / 'summary.md'
    result = subprocess.run(cmd + ['summary', '--projects', str(project), '--output', str(out)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert '## Tissue mixtures' in out.read_text() and 'linearity_error: unitas' in out.read_text()


# Depth and imbalance on a mixture (#488): its contrasts between fractions against the expected changes

def evaluated_mixture(tmp_path):
    """A prepared, quantified and evaluated mixture project (`declaration`, both methods counting as the model gives,
    content ratio 0.5), with the cutadapt route's reports (100 reads sequenced and 87 kept per sample)."""
    ext, path, rows = declaration(tmp_path)
    project = tmp_path / 'project'
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]
    result = subprocess.run(cmd + ['prepare', '--dataset', str(path), '--project', str(project), '--ext-dir', str(ext)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    names = list(yaml.safe_load((project / 'config.yaml').read_text())['samples'])
    quant = project / 'data/tmp/smallrna/quant'
    counts = model_counts(0.5, rows)
    write(quant / 'unitas/mir_counts.tsv', '\t' + '\t'.join(names) + '\n' + ''.join(
        '{}\t{}\n'.format(NAMES[m], '\t'.join('{:g}'.format(counts[s][m]) for s in names)) for m in sorted(PROFILE)))
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,' + ','.join(names) + '\n' + ''.join(
        'mmu-{},{}\n'.format(NAMES[m], ','.join('{:g}'.format(counts[s][m]) for s in names)) for m in sorted(PROFILE)))
    write_records(quant, ['unitas', 'mirge3'], dump=unitas_dump(ext, 'mus_musculus', 'Mouse', {}))
    for name in names:
        write_cutadapt_reports(project / 'data/tmp/smallrna/filter/cutadapt_stages', name, 100)
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    return project


def mixture_outputs(base, high, low, found=True):
    """de.smk's outputs for the mixture set in `base` (its inputs written): every method fits the changes the model gives
    between fractions `high` and `low`, called when `found`, else none called; `rare` untested."""
    w = lambda f: mixture.weight(f, 0.5)
    accession = {m: 'MIMAT00000{}'.format(10 + i) for i, m in enumerate(sorted(PROFILE))}
    change = {'liver_only': math.log2((1 - w(high)) / (1 - w(low))) if high < 1 else -8.0,
              'spleen_only': math.log2(w(high) / w(low)) if low > 0 else 8.0, 'both': 0.0, 'more_spleen': 1.2, 'more_liver': -0.3}
    lines = ['\t'.join(de.RESULT_HEADER)]
    for row in ('unitas', 'mirge3'):
        for test in de.TESTS:
            for m in sorted(PROFILE):
                called = found and abs(change.get(m, 0)) > 0.25
                fit = ('NA', 'NA', 'NA') if m == 'rare' else (change[m], 0.001, 0.01) if called else (change[m], 0.5, 0.9)
                lines.append('\t'.join([row, test, accession[m]] + [str(v) for v in fit]))
    text = '\n'.join(lines) + '\n'
    write(base / 'results.tsv', text)
    record = {'md5_' + Path(name).stem: de.md5((base / 'input' / name).read_text()) for name in de.INPUTS}
    record.update(md5_script=de.md5(de.SCRIPT.read_bytes()), md5_results=de.md5(text), R='R version 4.5.3', edgeR='4.10.7', limma='3.68.5',
                  formula='~subject + group', min_mean='10', edger='exact test', edger_glm='GLM test', voom='voom')
    write(base / 'record.tsv', 'field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))


def test_the_fixed_arm_shift_is_the_mirna_share_of_filtered_reads_between_fractions():
    """#488: the `fixed` arm divides by filtered reads; between two fractions that adds the log2 ratio of the miRNA share
    of filtered reads, mean over the subjects, to every fold change."""
    design, counts, _ = scored()
    spec = {'samples': samples()}
    between = [(s, numerator, denominator) for s, numerator, denominator, inner in benchmark.mixture_de_sets(spec, design) if inner]
    kept = {sample: 200000.0 if fraction >= 0.5 else 100000.0 for sample, fraction in
            ((row['sample'], row['fraction']) for row in samples())}   # twice the filtered reads from 50% on
    shift = benchmark.mixture_fixed_shift({'m': counts}, design, between, kept, ['m'])['m']
    assert shift['.75_vs_25'] == pytest.approx(-1.0, abs=0.001) and shift['.50_vs_25'] == pytest.approx(-1.0, abs=0.001)
    assert shift['.75_vs_50'] == pytest.approx(0.0, abs=0.001)   # the same miRNA reads and filtered reads at both


FRACTION = {'all': (1.0, 0.0), 'contrast.50_vs_25': (0.5, 0.25), 'contrast.75_vs_25': (0.75, 0.25), 'contrast.75_vs_50': (0.75, 0.5)}


def test_depth_de_scores_a_mixture_against_its_expected_changes(tmp_path):
    """`depth --de` on a mixture (#488): the mixture's sets at each depth, the full depth's inputs the evaluation's; the
    contrasts between fractions scored against each route's expected changes from the unthinned counts."""
    project = evaluated_mixture(tmp_path)
    bench, out = project / 'benchmark', tmp_path / 'depth'
    command = [sys.executable, str(BENCH / 'benchmark.py'), 'depth', '--project', str(project), '--depths', '20', '--output', str(out)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert (out / 'de/levels.tsv').read_text() == de.levels_text([(d, name) for d in ('full', '20') for name in FRACTION])
    for name in de.INPUTS:
        assert (out / 'de/full/all/input' / name).read_text() == (bench / 'de/input' / name).read_text()
        assert (out / 'de/full/contrast.75_vs_25/input' / name).read_text() == (bench / 'de/contrast/75_vs_25/input' / name).read_text()
    for depth in ('full', '20'):
        for name, (high, low) in FRACTION.items():
            mixture_outputs(out / 'de' / depth / name, high, low, found=depth == 'full')
    result = subprocess.run(command + ['--de'], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    tsv = (out / 'depth_de.tsv').read_text().splitlines()
    assert 'unitas\tfull\tbetween\tedger_glm_called\t1.0' in tsv and 'unitas\t20\tbetween\tedger_glm_called\t0.0' in tsv
    assert 'mirge3\tfull\tbetween\tvoom_wrong_direction\t0.0' in tsv and any(line.startswith('unitas\t20\tbetween\tscored\t') for line in tsv)
    assert 'unitas\t20\tall\tedger_recall\t0.0' in tsv   # the pure samples against the full depth's calls
    report = (out / 'depth.md').read_text()
    assert '## Between the mixtures, against the expected changes' in report and '| full | 1.00 (1.00) / 0.00 / ' in report
    assert '**The pure samples (2 100 and 2 0 samples of 2 subjects)**' in report and '**scored**' in report


def test_imbalance_scores_a_mixture_against_its_expected_changes(tmp_path):
    """`imbalance` on a mixture (#488): `tied` thins each set's lower fraction; the contrasts between fractions in each
    design and arm against the expected changes."""
    project = evaluated_mixture(tmp_path)
    bench, out = project / 'benchmark', tmp_path / 'imbalance'
    command = [sys.executable, str(BENCH / 'benchmark.py'), 'imbalance', '--project', str(project), '--ratios', '2', '--output', str(out)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    folders = ['{}_{}/{}'.format(scenario, arm, name) for scenario in ('none', 'independent-2', 'tied-2') for name in FRACTION
               for arm in ('workflow', 'fixed', 'tmm')]
    assert (out / 'de/levels.tsv').read_text() == de.levels_text([folder.split('/') for folder in folders])
    for name in de.INPUTS:
        assert (out / 'de/none_workflow/contrast.75_vs_25/input' / name).read_text() == (bench / 'de/contrast/75_vs_25/input' / name).read_text()
    sizes = lambda folder: [float(v) for v in (out / 'de' / folder / 'input/library_sizes.tsv').read_text().splitlines()[1].split('\t')[1:]]
    tied = [t / f for t, f in zip(sizes('tied-2_workflow/contrast.75_vs_25'), sizes('none_workflow/contrast.75_vs_25'))]
    assert tied[:2] == [1.0, 1.0] and all(0.45 < r < 0.55 for r in tied[2:])   # A_075, B_075 kept; A_025, B_025 halved
    assert sizes('tied-2_fixed/contrast.75_vs_25') == [87.0, 87.0, 43.5, 43.5]
    assert 'tied-2\tcontrast.75_vs_25\t87.0\t43.5\t2.0' in (out / 'imbalance_depths.tsv').read_text()
    for folder in folders:
        scenario, name = folder.split('_', 1)[0], folder.split('/')[1]
        mixture_outputs(out / 'de' / folder, *FRACTION[name], found=scenario == 'none')
    result = subprocess.run(command + ['--de'], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    tsv = (out / 'imbalance.tsv').read_text().splitlines()
    assert 'unitas\tnone\tworkflow\tbetween\tedger_called\t1.0' in tsv and 'unitas\ttied-2\tfixed\tbetween\tedger_called\t0.0' in tsv
    assert 'mirge3\ttied-2\ttmm\tall\tedger_recall\t0.0' in tsv
    report = (out / 'imbalance.md').read_text()
    assert '## Between the mixtures, against the expected changes' in report and '| tied-2 | tmm | 0.00 (0.00) / - / ' in report
    assert 'it is scored against the expected changes shifted by that ratio' in report and '| 75 against 25 | +0.00 | +0.00 |' in report
