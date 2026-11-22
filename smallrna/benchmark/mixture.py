"""Tissue mixtures (#451): biological samples of known structure. Total RNA of two tissues (`mixture.components`, e.g.
liver and spleen) of one subject is mixed at known fractions; the composition of no sample is known, but each mixture
is the weighted sum of its subject's two pure samples. So every method gets a truth from its own counts: its counts of
a mixture are predicted from its counts of the subject's pure samples, and the residual scores how linear and
consistent its counting is across samples. It does not score assignment (a method that misassigns the same way in
every sample still passes); the pools and the synthetic dataset do that.

The model. `fraction` is the share of the second component in the total RNA. The tissues hold different amounts of
miRNA per unit of RNA, and the library preparation recovers miRNAs unequally, so the share of a mixture's miRNA reads
that comes from the second component is not `fraction` but

    w(fraction) = fraction * k / (fraction * k + 1 - fraction)

with a content ratio k (the second component's miRNA reads per unit of RNA over the first's, as a method counts them)
for each subject: its two tissue samples have their own RNA content, and a per-miRNA bias acts on each subject's own
profiles. A per-miRNA bias that is the same in every sample changes k and nothing else: the observed proportions of a
mixture stay w * (pure second) + (1 - w) * (pure first) of the observed pure proportions. k is fitted per method and
subject (content_ratio): each mixture's weight by least squares of its proportions on the subject's two pure profiles,
turned into k, and the median over the subject's mixtures. It is fitted on the mixtures that are then scored: one
number per subject against the residuals of all its miRNAs.

Samples: a subject's series is its two pure samples (fraction 0 and 1) and its mixtures, one sample per fraction (the
lowest `replicate`). Further replicates of a subject and fraction are technical replicates: they are not in the series,
and give the noise two libraries of one RNA have (replicate_error), the floor for the prediction's error. A subject
without both pure samples (a pool of several subjects' RNA) has no series.

Metrics per method (scores), over miRNA-sample pairs of the mixtures with at least MIN_READS predicted reads:
- linearity_error: median |log2(observed / predicted)| (half a read added to both);
- within_1_5fold: the share within 1.5-fold;
- monotonic: of the miRNA-subject series of tissue-specific miRNAs (at least SPECIFIC log2 apart between the pure
  samples, MIN_READS predicted in every mixture), the share whose proportions never move against the trend over the
  fractions;
- replicate_error: median |log2 ratio| of the proportions of technical replicates, miRNAs with at least MIN_READS reads
  on average in the pair;
- content_ratio and content_ratio_spread (the median of the subjects' k, and the median absolute log2 deviation of
  the subjects' k from it): a description of the fit, not a score.

Differential expression (expected_changes): between two fractions a method's expected log2 fold change of a miRNA is
the mean over the subjects of log2(prediction at one / prediction at the other). benchmark.py tests the contrasts
between the mixtures' fractions with de.R, where the expectation comes from samples outside the contrast (the pure
ones), and the pure samples against each other as a description.
"""
import collections
import math

import numpy as np

MIN_READS = 10            # predicted (or mean) reads of a miRNA-sample pair that is scored
SPECIFIC = 2.0            # log2 difference between the pure samples of a tissue-specific miRNA
WITHIN = math.log2(1.5)
MIN_EXPECTED = 0.2        # least |expected log2 fold change| of a miRNA scored in a contrast: below it, no direction to score
BINS = ((10, 100), (100, 1000), (1000, None))   # predicted reads, for the error by abundance
KEYS = ('linearity_error', 'within_1_5fold', 'monotonic', 'replicate_error', 'content_ratio', 'content_ratio_spread')
SCORE_KEYS = KEYS[:4]     # with a bootstrap interval


class Design:
    """The samples of a mixture declaration. `series`: {subject: [(fraction, sample)]} by fraction, the subjects with
    both pure samples; `replicates`: [(sample, its technical replicate)]; `unused`: samples outside both."""

    def __init__(self, samples):
        rows = []
        for s in samples:
            fraction, subject = s.get('fraction'), s.get('subject')
            if isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not 0 <= fraction <= 1:
                raise ValueError('sample {}: fraction must be a number from 0 to 1, got {!r}'.format(s.get('sample'), fraction))
            if subject is None or not str(subject).strip():
                raise ValueError('sample {}: needs a subject'.format(s.get('sample')))
            replicate = s.get('replicate', 1)
            if isinstance(replicate, bool) or not isinstance(replicate, int) or replicate < 1:
                raise ValueError('sample {}: replicate must be a whole number of 1 or more, got {!r}'.format(s.get('sample'), replicate))
            rows.append((str(subject), float(fraction), replicate, s['sample']))
        by_key = collections.defaultdict(list)
        for subject, fraction, replicate, sample in sorted(rows):
            by_key[(subject, fraction)].append((replicate, sample))
        for (subject, fraction), members in by_key.items():
            if len({replicate for replicate, _ in members}) != len(members):
                raise ValueError('subject {} has two samples at fraction {:g} with the same replicate number: {}'.format(
                    subject, fraction, ', '.join(sample for _, sample in members)))
        subjects = sorted({subject for subject, _ in by_key})
        self.series, self.replicates, used = {}, [], set()
        for subject in subjects:
            fractions = sorted(fraction for s, fraction in by_key if s == subject)
            if 0.0 in fractions and 1.0 in fractions and len(fractions) > 2:
                self.series[subject] = [(fraction, by_key[(subject, fraction)][0][1]) for fraction in fractions]
                used.update(sample for _, sample in self.series[subject])
        for (subject, fraction), members in by_key.items():
            self.replicates += [(members[0][1], sample) for _, sample in members[1:]]
            used.update(sample for _, sample in members[1:])
        self.unused = sorted(sample for _, _, _, sample in rows if sample not in used)
        if len(self.series) < 2:
            raise ValueError('a mixture dataset needs at least two subjects with both pure samples (fraction 0 and 1) and a mixture, '
                             'got {}'.format(len(self.series)))
        self.fractions = sorted({fraction for members in self.series.values() for fraction, _ in members})
        self.samples = [sample for _, _, _, sample in rows]

    def at(self, fraction):
        """[(subject, sample)] of the series' samples at a fraction."""
        return [(subject, sample) for subject, members in self.series.items() for f, sample in members if f == fraction]


def label(fraction):
    """A fraction as a group label and part of a folder name: its percentage."""
    return '{:g}'.format(100 * fraction)


def contrasts(design):
    """[(numerator fraction, denominator fraction, between mixtures)]: first the pure samples, second component against
    first; then every pair of the mixtures' fractions, the higher against the lower. A subject is in a contrast when
    its series has both fractions."""
    inner = [f for f in design.fractions if 0 < f < 1]
    return [(1.0, 0.0, False)] + [(b, a, True) for i, a in enumerate(inner) for b in inner[i + 1:]]


def weight(fraction, k):
    """The share of a mixture's miRNA reads from the second component, at content ratio k."""
    return fraction * k / (fraction * k + 1 - fraction)


class Counts:
    """A method's counts over the design's samples: `matrix` (miRNA x sample), each sample's total and proportions."""

    def __init__(self, counts, samples, features):
        self.column = {sample: j for j, sample in enumerate(samples)}
        self.features = list(features)
        self.matrix = np.array([[counts.get(sample, {}).get(feature, 0.0) for sample in samples] for feature in features], dtype=float)
        self.matrix = self.matrix.reshape(len(self.features), len(samples))
        self.total = self.matrix.sum(axis=0)
        self.proportion = np.divide(self.matrix, self.total, out=np.zeros_like(self.matrix), where=self.total > 0)

    def of(self, sample):
        return self.proportion[:, self.column[sample]]


def mixture_weight(mixed, first, second):
    """The weight of `second` in `mixed` (proportions), by least squares: mixed = w * second + (1 - w) * first."""
    difference = second - first
    scale = float(difference @ difference)
    return float(difference @ (mixed - first)) / scale if scale > 0 else float('nan')


def content_ratio(counts, design):
    """({subject: k}, {sample: its own ratio}) of a method: per subject the median of its mixtures' ratios
    w (1 - f) / (f (1 - w)). A mixture without reads, or whose weight is not between 0 and 1, has no ratio, and a subject
    without such a mixture, or without reads in a pure sample, no k."""
    ks, ratios = {}, {}
    for subject, members in design.series.items():
        if min(counts.total[counts.column[sample]] for _, sample in members) <= 0:
            continue
        first, second = counts.of(members[0][1]), counts.of(members[-1][1])
        own = []
        for fraction, sample in members[1:-1]:
            w = mixture_weight(counts.of(sample), first, second)
            if 0 < w < 1:
                ratios[sample] = w * (1 - fraction) / (fraction * (1 - w))
                own.append(ratios[sample])
        if own:
            ks[subject] = float(np.median(own))
    return ks, ratios


def predicted(counts, members, fraction, k):
    """The predicted proportions of a subject's sample at `fraction`, from its pure samples (`members`: its series)."""
    w = weight(fraction, k)
    return w * counts.of(members[-1][1]) + (1 - w) * counts.of(members[0][1])


class Scores:
    """A method's residuals, kept as arrays so that a bootstrap resample is cheap (scores)."""

    def __init__(self, counts, design):
        self.ks, self.ratios = content_ratio(counts, design)
        self.k = float(np.median(list(self.ks.values()))) if self.ks else float('nan')
        self.spread = float(np.median([abs(math.log2(k / self.k)) for k in self.ks.values()])) if self.ks else float('nan')
        n = len(counts.features)
        self.subjects = sorted(design.series)
        columns, owner = [], []
        self.reads = []
        for subject in self.subjects:
            members = design.series[subject]
            for fraction, sample in members[1:-1]:
                expected = predicted(counts, members, fraction, self.ks[subject]) * counts.total[counts.column[sample]] \
                    if subject in self.ks else np.full(n, np.nan)
                observed = counts.matrix[:, counts.column[sample]]
                with np.errstate(invalid='ignore'):
                    residual = np.where(expected >= MIN_READS, np.log2((observed + 0.5) / (expected + 0.5)), np.nan)
                columns.append(residual)
                self.reads.append(expected)
                owner.append(subject)
        self.residual = np.column_stack(columns) if columns else np.empty((n, 0))   # miRNA x mixture sample
        self.reads = np.column_stack(self.reads) if columns else np.empty((n, 0))
        self.owner = owner
        trend = []
        for subject in self.subjects:
            members = design.series[subject]
            series = np.column_stack([counts.of(sample) for _, sample in members])
            totals = np.array([counts.total[counts.column[sample]] for _, sample in members])
            first, second = series[:, 0], series[:, -1]
            floor = 0.5 / max(float(totals.min()), 1.0)
            direction = np.log2((second + floor) / (first + floor))
            mixed = np.column_stack([predicted(counts, members, fraction, self.ks[subject]) * counts.total[counts.column[sample]]
                                     for fraction, sample in members[1:-1]]) if subject in self.ks else np.zeros((n, 1))
            eligible = (np.abs(direction) >= SPECIFIC) & (mixed >= MIN_READS).all(axis=1)
            steps = np.diff(series, axis=1) * np.sign(direction)[:, None]
            trend.append(np.where(eligible, (steps >= 0).all(axis=1).astype(float), np.nan))
        self.trend = np.column_stack(trend)   # miRNA x subject: 1 monotonic, 0 not, nan not tissue-specific or too few reads
        pairs = []
        for a, b in design.replicates:
            ca, cb = counts.matrix[:, counts.column[a]], counts.matrix[:, counts.column[b]]
            ta, tb = counts.total[counts.column[a]], counts.total[counts.column[b]]
            with np.errstate(invalid='ignore', divide='ignore'):
                ratio = np.log2((ca + 0.5) / ta) - np.log2((cb + 0.5) / tb) if ta > 0 and tb > 0 else np.full(n, np.nan)
            pairs.append(np.where((ca + cb) / 2 >= MIN_READS, ratio, np.nan))
        self.replicate = np.column_stack(pairs) if pairs else np.empty((n, 0))   # miRNA x replicate pair

    def values(self, features=None, subjects=None):
        """{metric: value} over the miRNAs (row indices; a resample repeats some) and subjects given; all by default."""
        rows = np.arange(self.residual.shape[0]) if features is None else np.asarray(features, dtype=int)
        subjects = self.subjects if subjects is None else subjects
        columns = [j for subject in subjects for j, owner in enumerate(self.owner) if owner == subject]
        residual = np.abs(self.residual[np.ix_(rows, columns)]) if len(rows) and columns else np.empty(0)
        residual = residual[~np.isnan(residual)]
        trend = self.trend[np.ix_(rows, [self.subjects.index(s) for s in subjects])] if len(rows) else np.empty(0)
        trend = trend[~np.isnan(trend)]
        replicate = np.abs(self.replicate[rows]) if len(rows) else np.empty(0)
        replicate = replicate[~np.isnan(replicate)]
        tidy = lambda value: round(float(value), 4)
        return {'linearity_error': tidy(np.median(residual)) if residual.size else float('nan'),
                'within_1_5fold': tidy((residual <= WITHIN).mean()) if residual.size else float('nan'),
                'monotonic': tidy(trend.mean()) if trend.size else float('nan'),
                'replicate_error': tidy(np.median(replicate)) if replicate.size else float('nan'),
                'content_ratio': tidy(self.k) if not math.isnan(self.k) else float('nan'),
                'content_ratio_spread': tidy(self.spread) if not math.isnan(self.spread) else float('nan')}

    def by_abundance(self):
        """[(bin label, pairs, median |log2 residual|)] of the scored miRNA-sample pairs by predicted reads (BINS)."""
        out = []
        for low, high in BINS:
            inside = (self.reads >= low) & ((self.reads < high) if high else True) & ~np.isnan(self.residual)
            values = np.abs(self.residual[inside])
            out.append(('{}-{}'.format(low, high - 1) if high else '>= {}'.format(low), int(values.size),
                        round(float(np.median(values)), 4) if values.size else float('nan')))
        return out


def expected_changes(counts, design, ks, numerator, denominator, offsets=None):
    """{feature: expected log2 fold change} of a method between two fractions: the mean over the subjects with both
    fractions and a content ratio (`ks`: {subject: k}) of log2(prediction at numerator / prediction at denominator), for
    the miRNAs a subject's pure samples hold. Empty without such a subject. `offsets` ({subject: log2 shift}) is added to
    each subject's log ratios before the mean (another scale, e.g. the `fixed` arm of #488)."""
    logs = []
    for subject, members in design.series.items():
        fractions = {fraction for fraction, _ in members}
        if subject not in ks or numerator not in fractions or denominator not in fractions:
            continue
        up, down = predicted(counts, members, numerator, ks[subject]), predicted(counts, members, denominator, ks[subject])
        with np.errstate(divide='ignore', invalid='ignore'):
            logs.append(np.where((up > 0) & (down > 0), np.log2(up / down), np.nan) + (offsets or {}).get(subject, 0.0))
    if not logs:
        return {}
    logs = np.column_stack(logs)
    defined = (~np.isnan(logs)).sum(axis=1)
    with np.errstate(invalid='ignore'):
        mean = np.nansum(logs, axis=1) / np.maximum(defined, 1)
    return {feature: float(mean[i]) for i, feature in enumerate(counts.features) if 2 * defined[i] >= logs.shape[1]}
