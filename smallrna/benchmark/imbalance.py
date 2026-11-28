"""Differential expression under controlled depth imbalance (#488): a stress test of the normalisation arms.

Each DE set of an evaluated biological project (the contrast, each level of `within`, the null comparisons of #449) is
thinned unevenly from the evaluation's count tables, binomially as the depth step thins (depth.thin), and tested again
with de.R in normalisation arms (`--arms`, default `workflow,fixed,tmmwsp`):
- `workflow`: a sample's library size is every miRNA read the route reported in it, as the workflow scales counts;
- `fixed`: the sample's filtered reads (the cutadapt route's kept reads) times its thinning fraction, for every route;
- `tmm`, `rle`, `tmmwsp` (#527): edgeR's TMM, RLE (median of ratios) or TMMwsp factors on the workflow's library sizes.
`--ratios none` runs the `none` design alone: every arm on the unthinned sets, e.g. the null comparisons' false calls
per normalisation on real data (#527).

Designs, each at the given ratios r:
- `independent-<r>`: in each of the set's two groups, half the samples (a seeded draw) thinned to 1/r, the others
  kept, so the imbalance is unrelated to the group;
- `tied-<r>`: the set's second group (the contrast's denominator, a null split's `y`) thinned to 1/r, so the first
  group is r times deeper than it was, as the adenomyosis samples are in the GSE317995 plasma series; only the second
  group is thinned, so a contrast is tested with its imbalance in one direction;
- `none`: the set unthinned, the control.

The ratios act on the project's own depths, which may already differ between the groups (a proportional read cap
keeps them): the realised depths, the groups' median filtered reads times their fractions, are reported per design.

Thinning keeps each sample's composition and every arm's library sizes follow it, so no arm's expected fold changes
move: the designs show what the imbalance costs (power, and false calls from the shallower samples' extra variance
under each test), not a normalisation bias from composition.

Scores, per design, arm, route and test: for the null comparisons `null_called`, `null_any` and `null_p05` (de.py)
and `offset`, the median over the splits of the median log2 fold change of the tested miRNAs (every change is null,
so a shift is the normalisation's); for the contrast's sets, against the same arm's unthinned fits, `calls`, `recall`
and `precision` (depth.de_scores) and `fc_shift`, the median change of the log2 fold changes. A mixture (#451): its
pure samples as a contrast, and its contrasts between fractions against each route's expected changes, from the
unthinned counts (`called`, `called_shared`, `wrong_direction`, `fc_bias`, `unchanged_called`, `scored`), the best
real test of what an arm does, since there every miRNA has an expected direction and size. The expected changes are
on the scale of a route's miRNA reads, the `workflow` arm's; the `fixed` arm gets them shifted, per subject, by the log2 ratio of the
route's miRNA share of filtered reads between the two fractions (benchmark.mixture_fixed_shift); `tmm` is scored on
the workflow's scale, as TMM meeting a design where most miRNAs change. `tied` thins each contrast's lower fraction,
the samples with less of the second tissue.
"""
import math
import statistics
import zlib

import numpy as np

import de
import depth

ALL_ARMS = ('workflow',) + de.ARMS  # every arm `--arms` may name (#453, #527)
ARMS = 'workflow,fixed,tmmwsp'  # the default, with #527's preset (tmmwsp): every arm multiplies the scenario fits
NONE = 'none'
DESIGNS = ('independent', 'tied')
RATIOS = '1.5,4'
CONTRAST_KEYS = depth.DE_KEYS + ('fc_shift',)
NULL_KEYS = de.NULL_KEYS + ('offset',)


def parse_arms(text):
    """Arms from 'workflow,fixed,tmm': names of ALL_ARMS, none twice, in ALL_ARMS' order."""
    arms = [part.strip() for part in str(text).split(',')]
    unknown = [a for a in arms if a not in ALL_ARMS]
    if unknown or not arms:
        raise ValueError('arms: {} not one of {}'.format(', '.join(unknown) or 'none given', ', '.join(ALL_ARMS)))
    if len(set(arms)) != len(arms):
        raise ValueError('arms: {} lists an arm twice'.format(text))
    return [a for a in ALL_ARMS if a in arms]


def parse_ratios(text):
    """Ratios from '1.5,4': numbers above 1, ascending, none twice; 'none': no thinning, the `none` design alone."""
    if str(text).strip() == NONE:
        return []
    ratios = []
    for part in str(text).split(','):
        try:
            value = float(part)
        except ValueError:
            raise ValueError('ratios: {!r} is not a number'.format(part.strip()))
        if not value > 1 or math.isinf(value):
            raise ValueError('ratios: {!r} must be a finite number above 1 (the deeper samples have r times the reads)'.format(part.strip()))
        ratios.append(value)
    if len(set(ratios)) != len(ratios):
        raise ValueError('ratios: {} lists a ratio twice'.format(text))
    return sorted(ratios)


def scenarios(ratios):
    """[(name, design, ratio)]: `none` first, then each design at each ratio."""
    return [(NONE, None, 1.0)] + [('{}-{:g}'.format(d, r), d, r) for r in ratios for d in DESIGNS]


def denominator(model):
    """The second group of a DE set's model.tsv text (de.model_text)."""
    fields = [line.split('\t', 1) for line in model.splitlines()[1:]]
    return next(value for field, value in fields if field == 'denominator')


def fractions(design, ratio, members, second, rng):
    """{sample: share of its reads kept} of a set's `members` ([(sample, group)]) under a design: `tied`, the
    `second` group's samples 1/ratio; `independent`, in each group the first half (rounded down) of its samples in a
    permutation drawn from `rng`, in samples.tsv order; None (`none`), every sample 1."""
    out = {sample: 1.0 for sample, _ in members}
    if design == 'tied':
        out.update({sample: 1 / ratio for sample, group in members if group == second})
    elif design == 'independent':
        for group in sorted({g for _, g in members}):
            samples = [sample for sample, g in members if g == group]
            out.update({samples[i]: 1 / ratio for i in rng.permutation(len(samples))[:len(samples) // 2]})
    elif design is not None:
        raise ValueError('unknown design {!r}'.format(design))
    return out


def generator(seed, scenario, folder):
    """The random generator of a set's draw under a scenario: seeded with `seed` and both names, so it does not depend
    on the order in which the sets are drawn."""
    return np.random.default_rng([seed, zlib.crc32('{}/{}'.format(scenario, folder).encode())])


def folder(scenario, arm, set_folder):
    """A set's folder under the DE root, `<scenario>_<arm>/<set>` (depth.de_set_name): one column and level of de.smk's
    levels.tsv."""
    return '{}_{}/{}'.format(scenario, arm, depth.de_set_name(set_folder))


def shifts(results, reference, features, rows):
    """{row: {<test>_fc_shift: median of log2fc - the reference's over the features both test}} (NaN without any)."""
    out = {}
    for row in rows:
        out[row] = {}
        for test in de.TESTS:
            pairs = [(results[(row, test, f)][0], reference[(row, test, f)][0]) for f in features]
            out[row][test + '_fc_shift'] = de.median([a - b for a, b in pairs if a is not None and b is not None])
    return out


def null_offsets(splits, rows):
    """{row: {<test>_offset: the median over `splits` ([(results, features)]) of the median log2 fold change of the
    tested features}} (NaN without any)."""
    out = {}
    for row in rows:
        out[row] = {}
        for test in de.TESTS:
            medians = [statistics.median(v) for v in ([results[(row, test, f)][0] for f in features if results[(row, test, f)][0] is not None]
                                                      for results, features in splits) if v]
            out[row][test + '_offset'] = de.median(medians)
    return out


def realised_rows(realised):
    """imbalance_depths.tsv's rows (design, set, first, second, ratio) from `realised` [(design, set folder, median
    filtered reads times fraction of the set's first group, of its second)]."""
    return [[scenario, depth.de_set_name(folder), round(first, 1), round(second, 1), round(first / second, 3) if second else '']
            for scenario, folder, first, second in realised]


def realised_lines(realised, md_table):
    """imbalance.md's table of realised depths: per design, the contrast's sets' medians and their ratio, and the null
    splits' ratio, median and range over the splits."""
    names = list(dict.fromkeys(folder for _, folder, _, _ in realised if not folder.startswith(de.NULL + '/')))
    has_nulls = any(folder.startswith(de.NULL + '/') for _, folder, _, _ in realised)
    table = []
    for scenario in dict.fromkeys(scenario for scenario, _, _, _ in realised):
        mine = [(folder, first, second) for sc, folder, first, second in realised if sc == scenario]
        cells = ['{:,.0f} / {:,.0f} ({:.2f})'.format(first, second, first / second) for folder, first, second in mine
                 if folder in names]
        if has_nulls:
            ratios = [first / second for folder, first, second in mine if folder.startswith(de.NULL + '/')]
            cells.append('{:.2f} ({:.2f}–{:.2f})'.format(statistics.median(ratios), min(ratios), max(ratios)))
        table.append([scenario] + cells)
    header = ['design'] + [depth.de_set_name(name) for name in names] + (['null splits'] if has_nulls else [])
    return ['', '## Realised depths', '',
            'Filtered reads times the thinning fraction, median per group: first / second (their ratio); for the null splits, '
            'the ratio\'s median and range over the splits. The ratios act on these depths, which may differ already in `none`.',
            ''] + md_table(header, table)


def rows(contrasts, nulls, methods, between=()):
    """imbalance.tsv's rows (method, scenario, arm, set, metric, value): `contrasts` [(set, title, [(scenario, arm,
    {row: values})])]; `nulls` [(scenario, arm, {row: values})], the set `null`; `between` likewise, a mixture's
    contrasts between fractions, the set `between`."""
    out = [[m, scenario, arm, name, key, value] for name, _, points in contrasts for scenario, arm, values in points
           for m in methods for key, value in values[m].items()]
    out += [[m, scenario, arm, de.NULL, key, value] for scenario, arm, values in nulls for m in methods for key, value in values[m].items()]
    return out + [[m, scenario, arm, depth.BETWEEN, key, value] for scenario, arm, values in between for m in methods
                  for key, value in values[m].items()]


signed = depth.signed


def between_notes(arms, shifts, methods, md_table):
    """imbalance.md's lines on the scales the selected `arms` are scored on between a mixture's fractions: `shifts`
    [(contrast, [the `fixed` arm's shift per method])] (benchmark.mixture_fixed_shift), shown when `fixed` is selected."""
    text = ['The expected changes are on the scale of a route\'s miRNA reads, the `workflow` arm\'s.']
    if 'fixed' in arms:
        text.append('The `fixed` arm divides by filtered reads instead, which adds each subject\'s log2 ratio of the route\'s '
                    'miRNA share of filtered reads, higher fraction over lower, to that subject\'s fold changes; it is scored against '
                    'expected changes with that ratio added to each subject\'s log ratios before their mean, from the unthinned data '
                    '(below, the mean over the subjects).')
    factors = [a for a in arms if a in de.NORMALISATIONS]
    if factors:
        text.append('{} {} scored on the workflow\'s scale: {} rescale towards most miRNAs not changing, which does not hold '
                    'between two fractions, so {} `fc_bias` and calls include that shift, as on the ratiometric pools.'.format(
                        ', '.join('`{}`'.format(a) for a in factors), 'is' if len(factors) == 1 else 'are',
                        'its factors' if len(factors) == 1 else 'their factors', 'its' if len(factors) == 1 else 'their'))
    table = md_table(['contrast'] + methods, [[name] + [signed(v) for v in values] for name, values in shifts]) if 'fixed' in arms else []
    return ['', ' '.join(text), ''] + table


ARM_TEXT = {'workflow': '`workflow` ({library})', 'fixed': '`fixed` (the sample\'s filtered reads times its thinning fraction)',
            'tmm': '`tmm` (edgeR\'s TMM factors on the workflow\'s)', 'rle': '`rle` (edgeR\'s RLE factors, median of ratios, on the workflow\'s)',
            'tmmwsp': '`tmmwsp` (edgeR\'s TMMwsp factors on the workflow\'s)'}


def lines(name, contrasts, nulls, methods, ratios, seed, versions, md_table, between=(), notes=(), library='every miRNA read the route reported',
          arms=None):
    """imbalance.md's lines; `library`: what the `workflow` arm's library size counts (a mixture's differ); `arms` the arms
    run (default ARMS)."""
    share = depth.share_cell
    arms = parse_arms(ARMS) if arms is None else arms
    named = [ARM_TEXT[a].format(library=library) for a in arms]
    out = ['# Depth imbalance: {}'.format(name), '',
           'Each DE set thinned unevenly from the evaluation\'s count tables, binomially (seed {}), and tested again with de.R in '
           'normalisation arms (#488, #527): {}. `independent-<r>`: in each '
           'group of the set, half the samples thinned to 1/r; `tied-<r>`: the set\'s second group only thinned to 1/r, so a '
           'contrast is tested with the imbalance in one direction; `none`: unthinned. Ratios: {}. Thinning keeps each sample\'s composition and every arm\'s library sizes '
           'follow it, so the designs show what the imbalance costs (power, false calls from the shallower samples\' variance), '
           'not a bias from composition. Each set is scored over the miRNAs of its evaluation; a row tests those its filter keeps in '
           'the thinned counts (#527\'s preset: filterByExpr).'.format(
               seed, ', '.join(named[:-1]) + (' and ' if len(named) > 1 else '') + named[-1], ', '.join('{:g}'.format(r) for r in ratios) or 'none')]
    if nulls:
        out += ['', '## Null comparisons (#449)', '',
                'The splits of each group, every call false: `null_called` (the median share of tested miRNAs called), the share '
                'of splits with a call in brackets, and `offset`, the median over the splits of the median log2 fold change, a shift '
                'that comes from the normalisation.']
        for test in de.TESTS:
            out += ['', '**{}**'.format(de.TEST_NAMES[test]), '']
            out += md_table(['design', 'arm'] + methods, [[scenario, arm] + ['{} ({}) {}'.format(
                share(values[m][test + '_null_called']), share(values[m][test + '_null_any']), signed(values[m][test + '_offset']))
                for m in methods] for scenario, arm, values in nulls])
    for _, title, points in contrasts:
        out += ['', '## {}'.format(title), '',
                'Against the same arm\'s unthinned calls (BH-adjusted p of at most {}, in the same direction): `recall / precision '
                '(calls)`, then `fc_shift`, the median change of the log2 fold changes; `none` shows the calls.'.format(de.FDR)]
        for test in de.TESTS:
            out += ['', '**{}**'.format(de.TEST_NAMES[test]), '']
            out += md_table(['design', 'arm'] + methods, [[scenario, arm] + [
                '{} / {} ({}) {}'.format(share(values[m][test + '_recall']), share(values[m][test + '_precision']), values[m][test + '_calls'],
                                         signed(values[m][test + '_fc_shift'])) if test + '_recall' in values[m]
                else str(values[m][test + '_calls']) for m in methods] for scenario, arm, values in points])
    if between:
        out += depth.between_lines([([scenario, arm], values) for scenario, arm, values in between], ['design', 'arm'], methods, md_table,
                                   'each design and arm') + list(notes)
    return out + ['', versions]
