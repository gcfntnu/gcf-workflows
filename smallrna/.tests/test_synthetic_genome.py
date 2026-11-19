"""Genome-wide synthetic dataset (#412): the generator (synthetic_genome.py), its metrics (genome_metrics.py) and the
benchmark's prepare and evaluate in genome mode, on a small fake miRBase with each hard case in it."""
import collections
import copy
import gzip
import math
from pathlib import Path
import random
import subprocess
import sys

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / 'smallrna/benchmark'
sys.path.insert(0, str(BENCH))
import adapters  # noqa: E402
import benchmark  # noqa: E402
import fill  # noqa: E402
import genome_metrics as gm  # noqa: E402
import synthetic_genome as sg  # noqa: E402

ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'
SPEC = {
    'mode': 'genome', 'seed': 2, 'mirnas': 16,
    'abundance': {'zipf': 1.0, 'tail': {'mirnas': 3, 'reads': [1, 5]}},
    'near_identical': {'ratio': 100, 'partner_absent': 0.0},
    'reads': {'mean': 3000, 'cv': 0.2},
    'groups': {'control': 3, 'case': 3},
    'dispersion': {'asymptotic': 0.01, 'extra': 0.1},
    'structure': {'group': {'fraction': 0.25, 'log2': [-2, 2], 'include': ['near_identical']},
                  'nuisance': {'sd': 0.3, 'group_correlation': 0.3}},
    'coregulate': {'distance': 10000, 'shared': 0.8},
    'isomirs': {'canonical': 0.8, '3p_trim_1': 0.06, '3p_trim_2': 0.02, '3p_templated_1': 0.04, '3p_templated_2': 0.01,
                '3p_add_A': 0.03, '3p_add_U': 0.02, '5p_shift': 0.02},
    'errors': {'rate': 0.001, 'first_base': 0.02},
    'qualities': {'correct': {37: 0.9, 23: 0.08, 12: 0.02}, 'error': {12: 0.6, 2: 0.25, 23: 0.15}},
    'random_reads': {'per_million': 1000, 'length': [18, 20]},
    'probes': {'replicates': 1, 'carrier': {'mirnas': 3, 'reads': 20}, 'probe': {'mirnas': 3, 'reads': 10}},
}


def fake_mirbase(folder, n=40, seed=3):
    """mature.fa, hairpin.fa and genes.gff3 of n human miRNAs (miR-T1 ...), with the hard cases: T1 and T2 share their
    mature sequence; T3 and T4 differ at one position; T5 is T6 without its last base; T7 is at two loci; T8 and T9 lie
    within 10 kb; T10's hairpin has an A after the mature (its added A is templated)."""
    rng = random.Random(seed)
    rand = lambda k: ''.join(rng.choice('ACGT') for _ in range(k))
    matures = {i: rand(22) for i in range(1, n + 1)}
    matures[2] = matures[1]
    matures[4] = matures[3][:10] + ('A' if matures[3][10] != 'A' else 'C') + matures[3][11:]
    matures[6] = rand(23)
    matures[5] = matures[6][:-1]
    mature_fa, hairpin_fa, gff = [], [], ['##gff-version 3']
    for i in range(1, n + 1):
        after = 'A' + rand(19) if i == 10 else ('C' + rand(19))
        hairpin = rand(10) + matures[i] + after + rand(22) + rand(10)
        loci = [(i * 100000, '+')] + ([(i * 100000 + 500000, '-')] if i == 7 else [])
        if i == 9:
            loci = [(8 * 100000 + 3000, '+')]
        mature_fa.append('>hsa-miR-T{0}-5p MIMAT9{0:06d} Homo sapiens miR-T{0}-5p\n{1}'.format(i, matures[i].replace('T', 'U')))
        for k, (start, strand) in enumerate(loci):
            hp = 'MI9{:06d}'.format(i * 10 + k)
            hairpin_fa.append('>hsa-mir-T{}{} {} Homo sapiens mir-T{} stem-loop\n{}'.format(i, '-{}'.format(k + 1) if len(loci) > 1 else '', hp, i,
                                                                                            hairpin.replace('T', 'U')))
            gff.append('chr1\t.\tmiRNA_primary_transcript\t{}\t{}\t.\t{}\t.\tID={};Alias={};Name=hsa-mir-T{}'.format(
                start, start + len(hairpin) - 1, strand, hp, hp, i))
            gff.append('chr1\t.\tmiRNA\t{}\t{}\t.\t{}\t.\tID=MIMAT9{:06d}{};Alias=MIMAT9{:06d};Name=hsa-miR-T{}-5p;Derives_from={}'.format(
                start + 10, start + 10 + len(matures[i]) - 1, strand, i, '_1' if k else '', i, i, hp))
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'mature.fa').write_text('\n'.join(mature_fa) + '\n>mmu-miR-X MIMAT0000001 Mus musculus miR-X\nACGUACGUACGUACGUACGUAC\n')
    (folder / 'hairpin.fa').write_text('\n'.join(hairpin_fa) + '\n')
    (folder / 'genes.gff3').write_text('\n'.join(gff) + '\n')
    return [folder / 'mature.fa', folder / 'hairpin.fa', folder / 'genes.gff3'], matures


@pytest.fixture
def mirbase(tmp_path):
    return fake_mirbase(tmp_path / 'mirbase')


def generate(tmp_path, files, spec=SPEC):
    out = tmp_path / 'fastq'
    truth, samples, mirnas = sg.generate(copy.deepcopy(spec), *files, 'hsa-', ADAPTER, 0, lambda s: out / '{}_R1.fastq.gz'.format(s))
    return truth, samples, mirnas, out


def tables(truth, samples, mirnas):
    return ([dict(zip(sg.TRUTH_COLUMNS, r)) for r in truth], [dict(zip(sg.SAMPLE_COLUMNS, map(str, r))) for r in samples],
            [dict(zip(sg.MIRNA_COLUMNS, map(str, r))) for r in mirnas])


def test_universe_collapses_identical_sequences_and_finds_relatives(mirbase):
    files, matures = mirbase
    u = sg.Universe(*files, 'hsa-')
    assert u.alias['MIMAT9000002'] == 'MIMAT9000001' and u.mirnas['MIMAT9000001']['accessions'] == ['MIMAT9000001', 'MIMAT9000002']
    assert u.near['MIMAT9000003'] == {'MIMAT9000004'} and u.mirnas['MIMAT9000007']['loci'] == 2
    assert u.mirnas['MIMAT9000008']['cluster'] == u.mirnas['MIMAT9000009']['cluster'] != u.mirnas['MIMAT9000010']['cluster']
    assert 'MIMAT0000001' not in u.alias                      # another species' mature
    m = u.mirnas['MIMAT9000010']
    assert sg.variant(m, '3p_add_A', 1) is None and sg.variant(m, '3p_templated_1', 1) == matures[10] + 'A'
    assert sg.variant(m, '5p_shift', 1) == matures[10][1:] and sg.variant(m, '5p_shift', -1) == m['upstream'][-1] + matures[10]


def test_isomirs_equal_to_another_mature_are_left_out(mirbase):
    files, matures = mirbase
    u = sg.Universe(*files, 'hsa-')
    reads = sg.Reads(u, ['MIMAT9000006', 'MIMAT9000010'], SPEC, np.random.default_rng(1))
    assert '3p_trim_1' not in reads.reads['MIMAT9000006']    # it would be T5's mature
    assert reads.reads['MIMAT9000006']['canonical'] == matures[6] and '3p_add_A' not in reads.reads['MIMAT9000010']
    assert reads.shares['MIMAT9000006'][sg.READ_CLASSES.index('3p_trim_1')] == 0
    assert math.isclose(float(reads.shares['MIMAT9000006'].sum()), 1)
    # a read another miRNA could give is left out too, also when that miRNA is not in the truth
    assert '3p_trim_2' not in reads.reads['MIMAT9000006']    # T5's 3p_trim_1
    assert '3p_trim_1' not in sg.Reads(u, ['MIMAT9000005'], SPEC, np.random.default_rng(1)).reads['MIMAT9000005']  # T6's 3p_trim_2
    assert all(sg.possible_reads(m) <= reads.forbidden for m in u.mirnas.values())   # nor do errors or random reads give one


def test_an_added_base_that_any_hairpin_has_next_is_templated():
    seq = 'ACGTACGTACGTACGTACGTAC'
    m = {'sequence': seq, 'flanks': [('GGT', 'CTTA'), ('GGC', 'ATTA')], 'upstream': 'GGT', 'downstream': 'CTTA'}
    assert sg.variant(m, '3p_add_A', 1) is None and sg.variant(m, '3p_add_U', 1) == seq + 'T'
    assert sg.variant(m, '3p_templated_1', 1) == seq + 'C' and sg.variant(m, '5p_shift', -1) == 'T' + seq   # the first hairpin
    assert {seq + 'C', seq + 'AT', 'T' + seq, 'C' + seq, seq + 'A', seq[:-2], seq[1:]} <= sg.possible_reads(m)
    for unknown in (('GGC', ''), None):      # the mature ends the second hairpin, or is not found in it
        m['flanks'][1] = unknown
        assert sg.variant(m, '3p_add_A', 1) is None and sg.variant(m, '3p_add_U', 1) is None


class FamilyUniverse:
    """What choose() reads of a Universe: annotated miRNAs (one locus each) and their near-identical relatives."""

    def __init__(self, relatives, singles):
        self.near = collections.defaultdict(set)
        for a, b in relatives:
            self.near[a].add(b)
            self.near[b].add(a)
        self.mirnas = {r: {'name': r, 'loci': 1} for r in set(self.near) | set(singles)}

    def annotated(self):
        return sorted(self.mirnas)


def test_each_near_identical_pair_is_a_family_of_its_own():
    """Families like miR-548: a chain a-b-c-d, a triangle e-f-g and a pair h-i. Whatever is drawn, a chosen miRNA's only
    chosen relative is its partner, an absent relative belongs to one family, and a pair is never cut."""
    u = FamilyUniverse([('a', 'b'), ('b', 'c'), ('c', 'd'), ('e', 'f'), ('f', 'g'), ('e', 'g'), ('h', 'i')],
                       ['s{}'.format(k) for k in range(20)])
    spec = {'near_identical': {'ratio': 100, 'partner_absent': 0.5}, 'abundance': {'zipf': 1.0}}
    for seed in range(40):
        for n in (12, 3):
            chosen, base, roles, absent = sg.choose(u, dict(spec, mirnas=n), np.random.default_rng(seed))
            assert len(chosen) == len(set(chosen)) <= n and not absent & set(chosen) and set(base) == set(chosen)
            partner = {r: roles.get(r, ('', ''))[0] for r in chosen}
            for r in chosen:
                assert u.near[r] & set(chosen) <= {partner[r]} and (not partner[r] or partner[r] in chosen)
                assert bool(u.near[r]) == (r in roles)        # miRNAs with relatives come only as pairs or alone
            for p in absent:
                assert len({frozenset({r, partner[r]} - {''}) for r in u.near[p] if r in chosen}) == 1


def test_generation_is_deterministic_and_every_read_has_its_truth(tmp_path, mirbase):
    files, _ = mirbase
    first = generate(tmp_path / 'a', files)
    second = generate(tmp_path / 'b', files)
    assert first[:3] == second[:3]
    truth, samples, mirnas, out = first
    expected = collections.defaultdict(collections.Counter)
    for sample, _, _, _, seq, count, _ in truth:
        expected[sample][seq] += count
    for sample, *_ in samples:
        with gzip.open(out / '{}_R1.fastq.gz'.format(sample), 'rt') as fh:
            lines = fh.read().splitlines()
        inserts = collections.Counter(line[:-len(ADAPTER)] for line in lines[1::4])
        assert all(line.endswith(ADAPTER) for line in lines[1::4]) and inserts == expected[sample]
        assert all(len(q) == len(s) for s, q in zip(lines[1::4], lines[3::4]))
    kinds = collections.Counter(s[1] for s in samples)
    assert kinds == {'mixed': 6, 'probe': 8, 'control': 1}


def test_near_identical_pairs_and_counts_follow_the_design(mirbase):
    """The count model alone (no reads written): pairs at 100:1, and with almost no dispersion and no nuisance each
    miRNA's log2 case/control ratio is its true fold change plus the composition shift common to all miRNAs."""
    files, _ = mirbase
    spec = copy.deepcopy(SPEC)
    spec['reads'] = {'mean': 2000000, 'cv': 0.0}
    spec['dispersion'] = {'asymptotic': 0.0001, 'extra': 0.0}
    spec['structure']['nuisance'] = {'sd': 0.0, 'group_correlation': 0.0}
    rng = np.random.default_rng(spec['seed'])
    u = sg.Universe(*files, 'hsa-')
    chosen, base, roles, absent = sg.choose(u, spec, rng)
    assert {roles['MIMAT9000003'][1], roles['MIMAT9000004'][1]} == {'high', 'low'} and not absent & set(chosen)
    pair = sorted([base['MIMAT9000003'], base['MIMAT9000004']])
    assert math.isclose(pair[1] / pair[0], 100, rel_tol=1e-9)
    scores, loadings = sg.structure(chosen, base, u, spec, spec['groups'], rng)
    assert [s[0] for s in scores] == ['control'] * 3 + ['case'] * 3 and sum(1 for f, _ in loadings.values() if f) == 4
    table, library, phi = sg.counts(chosen, base, loadings, scores, spec, rng)
    case, control = [j for j, s in enumerate(scores) if s[0] == 'case'], [j for j, s in enumerate(scores) if s[0] == 'control']
    observed = {r: math.log2(sum(table[(j, r)] for j in case) / sum(table[(j, r)] for j in control))
                for r in chosen if base[r] > 2000}
    shift = float(np.median([v for r, v in observed.items() if not loadings[r][0]]))
    assert len(observed) > 4 and all(abs(v - shift - loadings[r][0]) < 0.1 for r, v in observed.items())


def test_probe_samples_hold_a_carrier_and_one_read_class(tmp_path, mirbase):
    files, _ = mirbase
    truth, samples, mirnas, _ = generate(tmp_path, files)
    classes = collections.defaultdict(set)
    errors = collections.defaultdict(list)
    for sample, cls, _, _, _, _, n_errors in truth:
        classes[sample].add(cls)
        errors[sample].append(n_errors)
    assert classes['P_control'] == {'carrier'}
    for cls in sg.ISOMIR_CLASSES:
        assert classes['P_{}_1'.format(cls)] == {'carrier', cls}
    assert classes['P_errors_1'] == {'carrier', 'canonical'}
    rows = [t for t in truth if t[0] == 'P_errors_1' and t[1] == 'canonical']
    assert rows and all(t[6] >= 1 for t in rows)
    carriers = {r[0] for r in mirnas if 'carrier' in r[13]}
    truth_ids = {r[0] for r in mirnas if r[4] > 0}
    assert len(carriers) == 3 and not carriers & truth_ids


def test_substituted_bases_get_low_qualities(mirbase):
    files, _ = mirbase
    u = sg.Universe(*files, 'hsa-')
    reads = sg.Reads(u, ['MIMAT9000011'], SPEC, np.random.default_rng(5))
    for _ in range(50):
        q = reads.quality(30, [0, 7])
        assert len(q) == 30 and {q[0], q[7]} <= {chr(33 + v) for v in (12, 2, 23)}
        assert set(q) <= {chr(33 + v) for v in (37, 23, 12, 2)}
    # the first cycle's raised error rate is on the insert's first base, or with UMIs on the UMI
    plain, umi = (sg.Reads(u, ['MIMAT9000011'], SPEC, np.random.default_rng(5), umi=k).error_rates(20) for k in (0, 4))
    assert plain[0] == 0.02 and set(plain[1:]) == {0.001} and set(umi) == {0.001}


def test_the_design_needs_two_groups(mirbase):
    files, _ = mirbase
    u = sg.Universe(*files, 'hsa-')
    for groups in ({'a': 2, 'b': 2, 'c': 2}, {'a': 4}, {'a': 1, 'b': 4}):
        with pytest.raises(ValueError):
            sg.structure(['MIMAT9000011'], {'MIMAT9000011': 1.0}, u, SPEC, groups, np.random.default_rng(1))


def test_t_distribution_and_bh():
    assert abs(gm.t_two_sided(2.228139, 10) - 0.05) < 1e-4 and gm.t_two_sided(0.0, 5) == pytest.approx(1)
    assert abs(gm.t_two_sided(-3.0, 4) - 0.0399420) < 1e-6   # scipy.stats.t: 2 * t.sf(3, 4)
    assert list(np.round(gm.bh([0.01, 0.04, 0.03, 0.2]), 4)) == [0.04, 0.0533, 0.0533, 0.2]


def test_metrics_on_perfect_spilling_and_isomir_dropping_methods(tmp_path, mirbase):
    files, _ = mirbase
    truth, samples, mirnas, _ = generate(tmp_path, files)
    T, S, M = tables(truth, samples, mirnas)
    roles = {m['mature_id']: (m['near_role'], m['near_identical'].split(',')) for m in M}
    perfect, canonical, spill = (collections.defaultdict(float) for _ in range(3))
    for r in T:
        if r['mature_id']:
            perfect[(r['sample'], r['mature_id'])] += r['count']
            if r['class'] in ('canonical', 'carrier') and not r['errors']:
                canonical[(r['sample'], r['mature_id'])] += r['count']
    for (s, i), c in perfect.items():
        role, partners = roles.get(i, ('', []))
        moved = 0.2 * c if role == 'high' else 0
        spill[(s, i)] += c - moved
        if moved:
            spill[(s, partners[0])] += moved
    common = [(m, s, '', i, c) for m, table in (('perfect', perfect), ('spill', spill), ('canonical', canonical)) for (s, i), c in table.items()]
    common.append(('perfect', 'M1', 'unknown-name', '', 3.0))
    ev = gm.Evaluation(common, T, S, M, ['perfect', 'spill', 'canonical'])
    p = ev.mixed_metrics('perfect', ev.units, ev.pairs, ev.de_units, ev.dev_units)
    assert p['recovered'] == 1 and p['precision'] == 1 and p['exact'] == 1 and p['near_spill'] == 0 and p['misassigned'] == 0
    assert p['de_sensitivity'] == ev.mixed_metrics('truth', ev.units, ev.pairs, ev.de_units, ev.dev_units)['de_sensitivity']
    assert 0 < ev.unplaced_share('perfect') < 0.01
    s = ev.mixed_metrics('spill', ev.units, ev.pairs, ev.de_units, ev.dev_units)
    assert 0.1 < s['near_spill'] < 0.3 and s['precision'] < 1
    c = ev.mixed_metrics('canonical', ev.units, ev.pairs, ev.de_units, ev.dev_units)
    assert c['recovered'] < 0.9 and c['precision'] == 1
    probes = ev.probes('canonical')
    assert all(v[0] == 0 for k, v in probes.items() if k != 'errors') and probes['errors'][0] == 0
    assert all(v[0] == pytest.approx(1) for v in ev.probes('perfect').values())
    r2_group, r2_nuisance, coexpression = ev.structure('truth')
    assert 0 <= r2_group <= 1 and 0 <= r2_nuisance <= 1
    # a method with the same count in every sample for all but one abundant miRNA has no co-expression to compare
    top = max(ev.units, key=lambda i: ev.t[ev.index[i]].sum())
    flat = [('flat', s, '', i, perfect[(s, i)] if i == top else 100.0) for i in ev.units for s in ev.mixed]
    assert math.isnan(gm.Evaluation(flat, T, S, M, ['flat']).structure('flat')[2])


def test_spill_onto_an_absent_relative_of_a_de_mirna_is_a_false_call(tmp_path, mirbase):
    # each row tests the miRNAs passing the count filter in its own counts, as an analysis would: a method that puts a
    # share of a DE miRNA's reads on an absent relative tests that relative, finds the same fold change and calls it
    files, _ = mirbase
    truth, samples, mirnas, _ = generate(tmp_path, files)
    T, S, M = tables(truth, samples, mirnas)
    perfect = collections.defaultdict(float)
    for r in T:
        if r['mature_id']:
            perfect[(r['sample'], r['mature_id'])] += r['count']
    common = [('perfect', s, '', i, c) for (s, i), c in perfect.items()]
    ev0 = gm.Evaluation(common, T, S, M, ['perfect'])
    source = next(i for i, v in sorted(ev0.de['truth'].items()) if v[0] and i in ev0.true_de)
    shadow = 'MIMAT0000001'   # a miRNA absent from the truth
    common += [('shadow', s, '', i, c) for (s, i), c in perfect.items()]
    common += [('shadow', s, '', shadow, 0.5 * c) for (s, i), c in perfect.items() if i == source]
    ev = gm.Evaluation(common, T, S, M, ['perfect', 'shadow'])
    assert shadow in ev.de_units and shadow not in ev.de['truth'] and ev.de['shadow'][shadow][0]
    p, s = (ev.mixed_metrics(name, ev.units, ev.pairs, ev.de_units, ev.dev_units) for name in ('perfect', 'shadow'))
    assert s['de_tested'] == p['de_tested'] + 1 and s['de_fdr'] > p['de_fdr'] and s['de_sensitivity'] == p['de_sensitivity']
    calls = {c[0]: c for c in ev.false_calls('shadow')}
    assert set(calls) == {c[0] for c in ev.false_calls('perfect')} | {shadow} and calls[shadow][1] == 0 and calls[shadow][2] > 0


def test_loader_and_fill_for_genome_mode(tmp_path):
    base = {'name': 'syn_genome', 'kind': 'synthetic', 'organism': 'homo_sapiens', 'libprepkit': 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4',
            'adapter': ADAPTER, 'methods': ['tinyrna'], 'data': 'benchmark/smallrna/syn_genome', 'synthetic': dict(SPEC),
            'mirbase': {'release': '23.1', 'url': 'https://example.org/mature.fa', 'sha256': 'a' * 64}}
    path = tmp_path / 'ds.yaml'
    path.write_text(yaml.safe_dump(base))
    with pytest.raises(SystemExit):                            # no hairpin and gff
        benchmark.load_dataset(path)
    base['mirbase'].update(hairpin={'url': 'https://example.org/hairpin.fa', 'sha256': 'b' * 64},
                           gff={'url': 'https://example.org/hsa.gff3', 'sha256': 'c' * 64})
    path.write_text(yaml.safe_dump(base))
    spec = benchmark.load_dataset(path)
    assert [p.name for _, _, _, p in fill.targets(spec, tmp_path)] == ['mature.fa', 'hairpin.fa', 'genes.gff3']
    downloads = []

    def fetch(url, algorithm, digest, path, log):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('x')
        downloads.append(url)
    assert fill.fill(spec, tmp_path, log=lambda *_: None, fetch=fetch) == 3 and len(downloads) == 3
    assert fill.fill(spec, tmp_path, log=lambda *_: None, fetch=fetch) == 0      # a second fill skips every file


def test_prepare_and_evaluate_in_genome_mode(tmp_path):
    ext = tmp_path / 'ext'
    files, _ = fake_mirbase(ext / 'benchmark/smallrna/syn_genome/mirbase/23.1')
    spec = {'name': 'syn_genome', 'kind': 'synthetic', 'organism': 'homo_sapiens', 'libprepkit': 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4',
            'adapter': ADAPTER, 'umi': 0, 'methods': ['tinyrna'], 'data': 'benchmark/smallrna/syn_genome', 'synthetic': copy.deepcopy(SPEC),
            'mirbase': {'release': '23.1', 'url': 'u', 'sha256': 'a' * 64, 'hairpin': {'url': 'u', 'sha256': 'b' * 64},
                        'gff': {'url': 'u', 'sha256': 'c' * 64}}}
    dataset = tmp_path / 'ds.yaml'
    dataset.write_text(yaml.safe_dump(spec))
    project = tmp_path / 'project'
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]
    result = subprocess.run(cmd + ['prepare', '--dataset', str(dataset), '--project', str(project), '--ext-dir', str(ext)],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert '6 mixed samples' in result.stdout and '9 probe samples' in result.stdout
    config = yaml.safe_load((project / 'config.yaml').read_text())
    assert 'db' not in config and len(config['samples']) == 15  # the standard references; mixed + probe + control
    bench = project / 'benchmark'
    truth = benchmark.read_tsv(bench / 'truth.tsv')
    counts = collections.defaultdict(float)
    for row in truth:
        if row['mature_id']:
            counts[(row['mature_id'], row['sample'])] += float(row['count'])
    samples = list(config['samples'])
    ids = sorted({i for i, _ in counts})
    table = ['Feature ID\tClassifier\tFeature Name\t' + '\t'.join(samples)]
    table += ['{}\tmiRNA\tname-{}\t'.format(i, i) + '\t'.join('{:g}'.format(counts.get((i, s), 0)) for s in samples) for i in ids]
    quant = project / 'data/tmp/smallrna/quant/tinyrna'
    quant.mkdir(parents=True)
    (quant / 'tinyrna_counts.tsv').write_text('\n'.join(table) + '\n')
    # the record quant_reference writes (#407): the miRBase release, bound to this count table
    (quant / 'reference.tsv').write_text('field\tvalue\tsource\nmirna_database\tmirbase\tx\nmirna_release\t23.1\tx\n'
                                         'table_sha256\t{}\tx\n'.format(adapters.table_sha256('tinyrna', quant.parent)))
    result = subprocess.run(cmd + ['evaluate', '--project', str(project), '--bootstrap', '20'], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    report = (bench / 'report.md').read_text()
    assert '| recovered | 1.0 | 1.0 |' in report and '| precision | 1.0 | 1.0 |' in report and '| near_spill | 0.0 | 0.0 |' in report
    assert '## Recovery by stratum' in report and '## Probe samples' in report and '| 3p_add_A |' in report
    assert (bench / 'metrics_ci.tsv').is_file() and (bench / 'per_mirna.tsv').read_text().startswith('mature_name\tmature_id\ttruth\ttinyrna\n')
    out = tmp_path / 'summary.md'
    result = subprocess.run(cmd + ['summary', '--projects', str(project), '--output', str(out)], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    row = next(line for line in out.read_text().splitlines() if line.startswith('| [project]'))
    assert '| 6 mixed, 9 probe |' in row and '| tinyrna |' in row
