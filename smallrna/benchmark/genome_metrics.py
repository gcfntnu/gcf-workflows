"""Metrics of the genome-wide synthetic dataset (#412), on the common table with identical sequences collapsed
(biological.collapse) and the truth written by synthetic_genome.py.

Mixed samples (`kind` mixed in samples.tsv), per method, over units: every miRNA in the truth or reported by any method.
- recovered: sum of min(truth, observed) / truth; precision: the same / observed, so reads on absent miRNAs and
  overcounting lower it; exact: share of truth miRNA-sample pairs counted exactly; misassigned: share of observed reads
  on miRNAs absent from the sample's truth.
- near_spill: over the near-identical pairs, reads above the truth on the low partner and any reads on an absent
  partner, relative to the high member's truth.
- deviation_correlation: miRNAs with at least DETECT truth reads on average; each miRNA's log2(count + 1) minus its mean
  over the samples, observed against truth, as one Pearson correlation.
- DE with de.R (#447, benchmark.py and de.py): the benchmark's DE step on these counts, scored over de_features; the
  metrics `<test>_*`. The true DE miRNAs are those designed DE with at least DE_MIN truth reads on average. (The
  built-in moderated t-test check that needed no container, `de_*`, was retired for #527's preset.)
- structure: the top two principal components of the centred log2 normalised counts (miRNAs with at least DETECT
  reads); r2 of each true factor (group, nuisance) regressed on them. coexpression: Spearman correlation between the
  method's and the truth's miRNA-miRNA correlations (miRNAs with at least COEXP reads).
- Strata: recovered among unique miRNAs, near-identical pairs, multi-locus matures, identical-sequence groups, and
  low-abundance miRNAs (under LOW reads per million).
Probe samples: per read class, the probe miRNAs' observed reads (minus the control's) / their truth.
"""
import collections

import numpy as np

DETECT = 20
DE_MIN = 10
COEXP = 50
LOW = 10
EXACT_TOLERANCE = 1e-6








def size_factors(matrix):
    """Median-of-ratios size factors of a miRNA x sample count matrix (rows positive in every sample)."""
    rows = matrix[(matrix > 0).all(axis=1)]
    if not len(rows):
        return np.ones(matrix.shape[1])
    log = np.log(rows)
    return np.exp(np.median(log - log.mean(axis=1, keepdims=True), axis=0))




def spearman(x, y):
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1]) if len(x) > 2 else float('nan')


class Evaluation:
    """Counts of the truth and each method over the mixed and probe samples, and what the metrics aggregate."""

    def __init__(self, common, truth, samples, mirnas, methods):
        self.methods = list(methods)
        self.mixed = [s['sample'] for s in samples if s['kind'] == 'mixed']
        self.samples = {s['sample']: s for s in samples}
        self.mirnas = {m['mature_id']: m for m in mirnas}
        self.truth = collections.defaultdict(float)
        self.decoy_reads = collections.defaultdict(float)  # decoy reads per sample (#412 PR 2): no miRNA, not random
        self.overlapped = collections.defaultdict(set)  # the matures an overlap probe's fragments overlap, per class
        for row in truth:
            if row['mature_id']:
                self.truth[(row['sample'], row['mature_id'])] += float(row['count'])
            elif row.get('source'):
                self.decoy_reads[row['sample']] += float(row['count'])
                if ':overlap:' in row.get('origin', ''):
                    self.overlapped[row['class']].add(row['origin'].rsplit(':', 1)[1])
        self.observed = collections.defaultdict(float)
        self.unplaced = collections.defaultdict(float)
        self.total = collections.defaultdict(float)
        for method, sample, _, mature_id, count in common:
            self.total[(method, sample)] += float(count)
            if mature_id:
                self.observed[(method, sample, mature_id)] += float(count)
            else:
                self.unplaced[(method, sample)] += float(count)
        truth_units = {i for (s, i) in self.truth if s in self.mixed}
        seen = {i for (m, s, i) in self.observed if s in self.mixed}
        self.units = sorted(truth_units | seen)
        self.t = np.array([[self.truth.get((s, i), 0.0) for s in self.mixed] for i in self.units])
        self.o = {m: np.array([[self.observed.get((m, s, i), 0.0) for s in self.mixed] for i in self.units]) for m in self.methods}
        self.index = {i: k for k, i in enumerate(self.units)}
        self.pairs = [(i, m['near_identical'].split(',')) for i, m in self.mirnas.items() if m['near_role'] in ('high', 'alone')]
        self.dev_units = [i for i in self.units if self.t[self.index[i]].mean() >= DETECT]
        self.true_de = {i for i in self.units if self.mirnas.get(i, {}).get('de') == '1' and self.t[self.index[i]].mean() >= DE_MIN}

    def counts(self, name):
        return self.t if name == 'truth' else self.o[name]

    def de_features(self):
        """[(unit, designed log2 fold change)] that de.py scores (#447): the true DE miRNAs with their fold change, and
        with 0 every other unit not designed DE (absent ones included). A miRNA designed DE with fewer than DE_MIN truth
        reads on average is left out: neither a change the design lets a row find nor an unchanged miRNA."""
        return [(i, float(self.mirnas[i]['log2fc']) if i in self.true_de else 0.0) for i in self.units
                if i in self.true_de or self.mirnas.get(i, {}).get('de') != '1']

    def mixed_metrics(self, name, units, pairs, dev_units):
        o = self.counts(name)
        rows = [self.index[i] for i in units]
        t, oo = self.t[rows], o[rows]
        positive = t > 0
        both = np.minimum(t, oo).sum()
        out = {'recovered': both / t.sum() if t.sum() else float('nan'),
               'precision': both / oo.sum() if oo.sum() else float('nan'),
               'exact': float((np.abs(oo - t)[positive] <= EXACT_TOLERANCE).mean()) if positive.any() else float('nan'),
               'misassigned': float(oo[~positive].sum() / oo.sum()) if oo.sum() else float('nan')}
        excess, high = 0.0, 0.0
        for i, partners in pairs:
            high += self.t[self.index[i]].sum() if i in self.index else 0.0
            for p in partners:
                if p in self.index:
                    k = self.index[p]
                    excess += np.maximum(o[k] - self.t[k], 0).sum()
        out['near_spill'] = excess / high if high else float('nan')
        dev = [self.index[i] for i in dev_units]
        if len(dev) > 2:
            lt, lo = np.log2(self.t[dev] + 1), np.log2(o[dev] + 1)
            dt, do = lt - lt.mean(axis=1, keepdims=True), lo - lo.mean(axis=1, keepdims=True)
            out['deviation_correlation'] = float(np.corrcoef(dt.ravel(), do.ravel())[0, 1])
        return out

    def strata(self, name, extra=None):
        """recovered per stratum; `extra` ({stratum: units}, e.g. MirGeneDB's, #493) adds strata, each shown even when
        empty."""
        groups = collections.defaultdict(list, {stratum: [] for stratum in (extra or {})})
        for i in self.units:
            m = self.mirnas.get(i)
            if not m or float(m['base_per_million']) <= 0:
                continue
            for stratum, units in (extra or {}).items():
                if i in units:
                    groups[stratum].append(i)
            if m['near_role']:
                groups['near_identical'].append(i)
            if int(m['loci']) > 1:
                groups['multi_locus'].append(i)
            if ',' in m['accessions']:
                groups['identical_sequences'].append(i)
            if float(m['base_per_million']) < LOW:
                groups['low_abundance'].append(i)
            if not m['near_role'] and int(m['loci']) == 1 and ',' not in m['accessions']:
                groups['unique'].append(i)
        o = self.counts(name)
        out = {}
        for stratum, units in groups.items():
            rows = [self.index[i] for i in units]
            t = self.t[rows].sum()
            out[stratum] = (float(np.minimum(self.t[rows], o[rows]).sum() / t) if t else float('nan'), len(units))
        return out

    def structure(self, name):
        """(r2 of the group factor, r2 of the nuisance factor) on the top two principal components; coexpression."""
        matrix = self.counts(name)
        rows = [self.index[i] for i in self.dev_units]
        if len(rows) < 3 or len(self.mixed) < 4:
            return float('nan'), float('nan'), float('nan')
        sf = size_factors(matrix[rows])
        y = np.log2(matrix[rows] / sf + 1)
        centred = y - y.mean(axis=1, keepdims=True)
        _, _, vt = np.linalg.svd(centred, full_matrices=False)
        pcs = np.column_stack([np.ones(len(self.mixed)), vt[0], vt[1]])
        r2 = []
        for key in ('factor_group', 'factor_nuisance'):
            v = np.array([float(self.samples[s][key]) for s in self.mixed])
            fit, *_ = np.linalg.lstsq(pcs, v, rcond=None)
            residual = v - pcs @ fit
            r2.append(float(1 - (residual ** 2).sum() / ((v - v.mean()) ** 2).sum()))
        co = [self.index[i] for i in self.units if self.t[self.index[i]].mean() >= COEXP]
        coexpression = float('nan')
        if len(co) > 3:
            obs = np.log2(matrix[co] + 1)
            keep = obs.std(axis=1) > 0
            if keep.sum() > 2:  # fewer varying miRNAs give too few correlations to compare
                truth_corr = np.corrcoef(np.log2(self.t[co] + 1))
                upper = np.triu_indices(int(keep.sum()), 1)
                coexpression = spearman(truth_corr[np.ix_(keep, keep)][upper], np.corrcoef(obs[keep])[upper])
        return r2[0], r2[1], coexpression

    def probes(self, name, members=None):
        """{class: (recovery, misplaced, probe miRNAs)}: observed reads on the class's probe miRNAs in its samples, minus
        the control's on them, over their truth; and the reads on other miRNAs (not carriers), minus the control's."""
        control = [s for s, row in self.samples.items() if row['kind'] == 'control']
        by_class = collections.defaultdict(list)
        for s, row in self.samples.items():
            if row['kind'] == 'probe':
                by_class[row['probe_class']].append(s)
        carriers = {i for i, m in self.mirnas.items() if 'carrier' in m['probe_classes']}
        out = {}
        for c, samples in sorted(by_class.items()):
            probe_ids = sorted({i for (s, i) in self.truth if s in samples} - carriers)
            chosen = probe_ids if members is None else members.get(c, probe_ids)
            obs = lambda s, i: self.observed.get((name, s, i), 0.0) if name != 'truth' else self.truth.get((s, i), 0.0)
            t = sum(self.truth.get((s, i), 0.0) for s in samples for i in chosen)
            got = sum(obs(s, i) for s in samples for i in chosen)
            background = sum(obs(cs, i) for cs in control for i in chosen) * len(samples) / max(len(control), 1)
            others = {i for (m, s, i) in self.observed if m == name and s in samples} - set(probe_ids) - carriers
            misplaced = sum(obs(s, i) for s in samples for i in others)
            misplaced -= sum(obs(cs, i) for cs in control for i in others) * len(samples) / max(len(control), 1)
            out[c] = ((got - background) / t if t else float('nan'), misplaced / t if t else float('nan'), len(probe_ids))
        return out

    def decoys(self, name, kind='decoy'):
        """{unit: (share, decoy reads, {miRNA: reads})} of the decoy probe samples (#412 PR 2; kind `decoy`, one class
        and source each; kind `overlap`, one class each): the method's reads on miRNAs other than the carrier, minus the
        control's on them (scaled to the number of samples), over the decoy reads they hold, and those reads per miRNA
        (above the control's)."""
        control = [s for s, row in self.samples.items() if row['kind'] == 'control']
        by_unit = collections.defaultdict(list)
        for s, row in self.samples.items():
            if row['kind'] == kind:
                by_unit[row['probe_class']].append(s)
        carriers = {i for i, m in self.mirnas.items() if 'carrier' in m['probe_classes']}
        obs = lambda s, i: self.observed.get((name, s, i), 0.0) if name != 'truth' else self.truth.get((s, i), 0.0)
        out = {}
        for unit, samples in sorted(by_unit.items()):
            seen = {i for (m, s, i) in self.observed if m == name and s in samples} - carriers if name != 'truth' else set()
            scale = len(samples) / max(len(control), 1)
            landed = {i: sum(obs(s, i) for s in samples) - sum(obs(cs, i) for cs in control) * scale for i in seen}
            reads = sum(self.decoy_reads[s] for s in samples)
            out[unit] = (sum(landed.values()) / reads if reads else float('nan'), reads, {i: v for i, v in landed.items() if v > 0})
        return out

    def overlaps(self, name):
        """{class: (on the overlapped miRNAs, on any miRNA, fragment reads, {miRNA: reads})} of the overlap probe samples
        (#412 PR 2): fragments of other-class molecules overlapping an embedded mature, counted as miRNAs (control
        subtracted) on the miRNAs they overlap, and on any; not false positives."""
        out = {}
        for cls, (share, reads, landed) in self.decoys(name, 'overlap').items():
            on = sum(v for i, v in landed.items() if i in self.overlapped[cls]) / reads if reads else float('nan')
            out[cls] = (on, share, reads, landed)
        return out

    def unplaced_share(self, name):
        total = sum(self.total[(name, s)] for s in self.mixed)
        return sum(self.unplaced[(name, s)] for s in self.mixed) / total if total else float('nan')
