"""Figures of the HTML summary (#338), each a PNG embedded in the page with a plain-language title and caption:

- how sure the per-route results are: a dot per route with its bootstrap interval (#413), for the equal-amount pools
  (within 2-fold of the exact-match count), one synthetic set (reads on the right miRNA) and the mixture (within
  1.5-fold of the mixing proportions);
- changes called where there are none: the null comparisons (#449) of a biological project, per test;
- sequencing depth: per abundance stratum the share of the full data's miRNAs still seen in fewer reads, and the new
  miRNAs per doubling of depth, from `benchmark depth` (#409) in a project's `benchmark/depth/`.

Projects are taken by what they declare and hold, as in the overview (overview.py); settings variants are left out.
A figure without a project is left out, and all of them without matplotlib.
"""
import base64
import io
import math

import overview

# a fixed colour per route, so a route keeps its colour across figures; other routes are grey
ROUTE_COLOR = {'unitas': '#2a78d6', 'mirge3': '#eb6834', 'tesmall': '#1baf7a', 'tinyrna': '#eda100', 'smrnaseq': '#e87ba4',
               'biocore_htseq': '#008300', 'biocore_miraligner': '#6250d6'}
OTHER = '#8a9199'
LABEL = {'unitas': 'UNITAS', 'mirge3': 'miRge3', 'tesmall': 'TEsmall', 'tinyrna': 'tiny-count', 'smrnaseq': 'nf-core smrnaseq',
         'biocore_htseq': 'biocore HTSeq', 'biocore_miraligner': 'biocore miraligner', 'mirflowz': 'mirflowz'}
INK, INK2, GRID = '#1d2125', '#5b6168', '#e3e6ea'
STRATUM_SHADES = ['#164a8a', '#2a78d6', '#7fb0ea', '#b9d3f3']  # more to fewer reads per million; the lowest stratum is orange
LOWEST = '#eb6834'
STYLE = {'font.family': 'DejaVu Sans', 'font.size': 9.5, 'axes.edgecolor': '#9aa1a8', 'axes.labelcolor': INK2, 'xtick.color': INK2,
         'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False}
NULL_TESTS = (('edger', 'edgeR exact test'), ('edger_glm', 'edgeR GLM'), ('voom', 'limma-voom'))


def available():
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


def encoded(fig):
    """The figure as base64 PNG; the figure is closed."""
    import matplotlib.pyplot as plt
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode('ascii')


def headline(ax, title, subtitle, lines=None):
    """A bold title over a grey subtitle, both left-aligned above the axes; `lines` the subtitle lines to make room for,
    so panels side by side align their titles."""
    lines = max(lines or 0, subtitle.count('\n') + 1)
    ax.annotate(subtitle, (0, 1), xycoords='axes fraction', xytext=(0, 6), textcoords='offset points', fontsize=8.2, color=INK2,
                va='bottom', linespacing=1.25)
    ax.annotate(title, (0, 1), xycoords='axes fraction', xytext=(0, 6 + 12.5 * lines + 5), textcoords='offset points', fontsize=10.5,
                color=INK, weight='bold', va='bottom')


def routes_of(projects):
    return list(dict.fromkeys(c for p in projects for c in p.columns + list(p.status) if c not in ('oracle', 'truth')))


def label(route):
    return LABEL.get(route, route)


def has_interval(project, metric):
    return any(m == metric for m, _ in project.intervals)


def interval_panels(projects):
    """[(project, metric, title, subtitle)] of the interval figure. The synthetic panel is one draw (intervals of
    different draws are not averaged): the lowest seed of the declaration with the most draws (one dataset name)."""
    panels = []
    pools = sorted((p for p in projects if p.kind == 'pool' and p.pool == 'equimolar' and has_interval(p, 'within_2fold')),
                   key=lambda p: (p.kit, p.path.name))
    for i, p in enumerate(pools):
        kit = overview.kit_label(p.kit)
        panels.append((p, 'within_2fold', 'Equal-amount pool, {}'.format(kit),
                       'Share of the pool\'s miRNAs whose count is within\n2-fold of the exact-match count of the same reads\n'
                       '(1 = every miRNA counted right)' if i == 0 else 'The same measure for the {}\nlibrary'.format(kit)))
    synthetic = [p for p in projects if p.kind == 'synthetic' and p.genome and p.value('de_sensitivity', 'truth')]
    shown = [p for p in synthetic if has_interval(p, 'recovered')]
    if shown:
        draws = {}
        for q in shown:
            draws.setdefault(q.spec['name'], []).append(q)
        name = min(draws, key=lambda n: (-len(draws[n]), n))
        p = min(draws[name], key=lambda q: (not isinstance(q.seed, int), q.seed if isinstance(q.seed, int) else 0, q.path.name))
        panels.append((p, 'recovered', 'Synthetic reads' + (' (seed {})'.format(p.seed) if p.seed is not None else ''),
                       'Share of the true miRNA reads the route\ncounts on the right miRNA (1 = all of them);\n'
                       'one of the {} synthetic set{} the overview averages'.format(len(synthetic), '' if len(synthetic) == 1 else 's')))
    mixtures = sorted((p for p in projects if p.kind == 'mixture' and has_interval(p, 'within_1_5fold')), key=lambda p: p.path.name)
    for i, p in enumerate(mixtures):
        panels.append((p, 'within_1_5fold', 'Mixture of known proportions' + (' ({})'.format(i + 1) if len(mixtures) > 1 else ''),
                       'Share of miRNAs whose level in each mixture is\nwithin 1.5-fold of what the known mixing\nproportions predict (1 = all of them)'))
    return panels


def interval_figure(projects, routes):
    panels = interval_panels(projects)
    if not panels:
        return None
    import matplotlib.pyplot as plt
    columns = min(4, len(panels))
    rows = math.ceil(len(panels) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(3.9 * columns, 4.4 * rows), sharey=True, squeeze=False)
    for ax, (project, metric, title, subtitle) in zip(axes.flat, panels):
        drawn = []
        for i, route in enumerate(routes):
            interval = project.intervals.get((metric, route))
            if interval is None:
                state = project.status.get(route) or ('no interval' if project.value(metric, route) is not None else '')
                ax.text(0.02, i, state, transform=ax.get_yaxis_transform(), va='center', fontsize=8, color=INK2, style='italic')
                continue
            value, low, high = interval
            color = ROUTE_COLOR.get(route, OTHER)
            ax.plot([low, high], [i, i], color=color, lw=2.4, solid_capstyle='round')
            ax.plot(value, i, 'o', color=color, ms=7.5, mec='white', mew=1.5, zorder=3)
            drawn.append(interval)
        headline(ax, title, subtitle, max(p[3].count('\n') + 1 for p in panels))
        ax.grid(axis='x', color=GRID)
        ax.set_xlabel('better →', color=INK2, loc='right', fontsize=8.5)
        if len(drawn) > 1 and max(lo for _, lo, _ in drawn) <= min(hi for _, _, hi in drawn):
            ax.text(0.5, -0.22, 'all intervals overlap', transform=ax.transAxes, ha='center',
                    fontsize=8.3, color=INK, style='italic')
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    axes[0][0].set_yticks(range(len(routes)), [label(r) for r in routes])
    axes[0][0].set_ylim(len(routes) - 0.5, -0.6)
    fig.tight_layout(h_pad=4)
    resamples = {(p.record.get('bootstrap') or {}).get('resamples') for p, *_ in panels}
    caption = ('Point: the measure; line: its 95% bootstrap interval{}, the range the value would likely take with another draw of '
               'the same samples. Routes whose lines do not overlap differ clearly. Data: {}. The other synthetic sets and the '
               'set-ratio pools are in the tables below.').format(
        ' ({:,} resamples of the samples)'.format(resamples.pop()) if len(resamples) == 1 and None not in resamples else '',
        ', '.join('`{}`{}'.format(p.spec['name'], ' seed {}'.format(p.seed) if p.kind == 'synthetic' and p.seed is not None else '')
                  for p, *_ in panels))
    return ('How sure are the results?', encoded(fig), caption)


def null_figure(projects, routes):
    nulls = [p for p in projects if p.kind == 'biological' and not p.lanes and p.has('edger_null_p05')]
    if not nulls:
        return None
    import matplotlib.pyplot as plt
    tests = [(k, t, name) for k, (t, name) in enumerate(NULL_TESTS) if any(p.has(t + '_null_p05') for p in nulls)]
    fig, axes = plt.subplots(1, len(nulls), figsize=(8.5 * len(nulls), 0.5 * len(routes) + 1.2), sharex=True, sharey=True, squeeze=False)
    for n, (ax, project) in enumerate(zip(axes[0], nulls)):
        for i, route in enumerate(routes):
            shown = False
            for k, test, _ in tests:  # a test keeps its marker and place when another is missing
                v = project.value(test + '_null_p05', route)
                if v is not None:
                    ax.plot(v, i + (k - 1) * 0.22, 'osD'[k], color=ROUTE_COLOR.get(route, OTHER), ms=6.5, mec='white', mew=1)
                    shown = True
            if not shown:
                ax.text(0.02, i, project.status.get(route, ''), transform=ax.get_yaxis_transform(), va='center', fontsize=8,
                        color=INK2, style='italic')
        ax.axvline(overview.NULL_LEVEL, color=INK, lw=1.2)
        ax.text(overview.NULL_LEVEL, -0.75, ' 5%: expected by chance', fontsize=8, color=INK)
        ax.set_xlim(left=0)
        ax.set_xlabel('share of tested miRNAs with p ≤ 0.05        better: at or left of the line ←', fontsize=8.5)
        splits = ((project.spec.get('de') or {}).get('null') or {}).get('splits')
        headline(ax, 'Changes called where there are none' + (' ({})'.format(n + 1) if len(nulls) > 1 else ''),
                 'The samples of each group split into two random halves{}: no real differences,\n'
                 'so about 5% of miRNAs pass p ≤ 0.05 by chance'.format(', {} times'.format(splits) if splits else ''))
        ax.grid(axis='x', color=GRID)
        handles = [plt.Line2D([], [], marker='osD'[k], ls='', color=INK2, label=name) for k, _, name in tests]
        ax.legend(handles=handles, fontsize=8, frameon=False, loc='lower right', title='test', title_fontsize=8)
    axes[0][0].set_yticks(range(len(routes)), [label(r) for r in routes])
    axes[0][0].set_ylim(len(routes) - 0.5, -1.0)
    caption = ('Each mark: the share of the tested miRNAs with p at most 0.05, pooled over all splits, for {} of differential '
               'expression. Right of the line: more changes called than chance gives. Data: {}.').format(
        'one test' if len(tests) == 1 else '{} common tests'.format(len(tests)),
        ', '.join('`{}`'.format(p.spec['name']) for p in nulls))
    return ('Changes called where there are none', encoded(fig), caption)


def depth_value(text):
    """Reads of a depth label (500k, 8M), or None for the full data."""
    if text == 'full':
        return None
    return float(text[:-1]) * {'k': 1e3, 'M': 1e6}[text[-1]] if text[-1] in 'kM' else float(text)


def stratum_label(name, prefix):
    low, _, high = name[len(prefix):].partition('-')
    return '≥{} CPM'.format(low[:-1]) if low.endswith('+') else '{}–{} CPM'.format(low, high)


def depth_figure(project):
    """The depth figure of one project from its benchmark/depth/depth.tsv, for the first declared route it has."""
    import matplotlib.pyplot as plt
    import depth as depth_series
    rows = overview.read_tsv(project.path / 'benchmark' / 'depth' / 'depth.tsv')
    methods = list(dict.fromkeys(r['method'] for r in rows))
    route = next((m for m in project.spec.get('methods', []) if m in methods), methods[0] if methods else None)
    if route is None:
        return None
    values = {(r['depth'], r['metric']): overview.number(r['value']) for r in rows if r['method'] == route}
    depths = sorted({d for d, _ in values if d != depth_series.FULL}, key=depth_value)
    order = depths + ([depth_series.FULL] if any(d == depth_series.FULL for d, _ in values) else [])
    if len(order) < 2:
        return None
    x = {d: depth_value(d) / 1e6 for d in depths}
    if depth_series.FULL in order:  # the samples' own depths differ: a labelled position past the last thinned depth
        x[depth_series.FULL] = x[depths[-1]] * 1.8
    prefix = depth_series.STRATUM_PREFIX
    strata = sorted({m for _, m in values if m.startswith(prefix)},
                    key=lambda m: -float(m[len(prefix):].partition('-')[0].rstrip('+')))
    gains = [m for m in ('gain_doubling_detected_1', 'gain_doubling_detected_10') if any(k == m for _, k in values)]
    if not strata or not gains:  # a series written before the strata and gains (#409): its depth.md has what it holds
        return None
    fig, axes = plt.subplots(1, 2, figsize=(13, 3.9))
    ax = axes[0]
    shades = STRATUM_SHADES[:max(0, len(strata) - 1)] + [LOWEST]
    for metric, color in zip(strata, shades):
        points = [(x[d], values.get((d, metric))) for d in order if values.get((d, metric)) is not None]
        ax.plot([p[0] for p in points], [p[1] for p in points], marker='o', ms=4, lw=2, color=color,
                label=stratum_label(metric, prefix))
    ax.set_xscale('log')
    ax.set_xticks([x[d] for d in order], ['{:g}'.format(x[d]) if d != depth_series.FULL else 'full' for d in order])
    ax.minorticks_off()
    ax.set_xlabel('million reads sequenced per sample')
    ax.set_ylabel('share still seen ({}+ reads)'.format(depth_series.SEEN))
    ax.set_ylim(0, 1.03)
    ax.grid(color=GRID)
    ax.axhline(0.9, color=INK2, ls=':', lw=1)
    ax.text(1.005, 0.9, '90%', transform=ax.get_yaxis_transform(), fontsize=7.5, color=INK2, va='center')
    if strata:
        ax.legend(frameon=False, fontsize=8.5, loc='lower right', title='abundance at full depth', title_fontsize=8.5)
    headline(ax, 'How deep to sequence to see a miRNA', 'MiRNAs found in the full data, grouped by abundance (reads per\n'
             'million, CPM): share still seen in fewer reads')
    ax = axes[1]
    steps = list(zip(order, order[1:]))
    for metric, color, offset, name in (('gain_doubling_detected_1', '#a9c8ee', -0.2, '1+ reads'),
                                        ('gain_doubling_detected_10', '#2a78d6', 0.2, '10+ reads')):
        heights = [values.get((high, metric)) for _, high in steps]
        bars = ax.bar([k + offset for k in range(len(steps))], [math.nan if h is None else h for h in heights], width=0.38,
                      color=color, label=name)
        if steps[-1][1] == depth_series.FULL:
            bars[-1].set_hatch('///')
            bars[-1].set_edgecolor('white')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(len(steps)), ['{}→{}'.format(a, b) for a, b in steps], rotation=35, ha='right', fontsize=8)
    ax.set_ylabel('new miRNAs per doubling of depth')
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(axis='y', color=GRID)
    headline(ax, 'What each doubling of depth adds', 'New miRNAs found each time the reads double, counted\nat 1+ and at 10+ reads; '
             'falling bars = saturating')
    fig.tight_layout()
    caption = ('Data: `{}` ({} samples), counted with {}; the other routes are in its depth.md. Fewer reads are simulated by '
               'dropping reads at random from the full data. Full: each sample\'s whole data; where the samples differ in depth, the '
               'last step (hatched) is approximate.').format(project.spec['name'], project.samples, label(route))
    return ('Sequencing depth', encoded(fig), caption)


def build(paths):
    """[(title, base64 PNG, caption)] of the projects at `paths`, [] when none applies; None without matplotlib."""
    if not available():
        return None
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    projects = [p for p in (overview.Project(path) for path in paths) if not p.variant]
    routes = routes_of(projects)
    out = []
    with plt.rc_context(STYLE):
        out += [f for f in (interval_figure(projects, routes), null_figure(projects, routes)) if f]
        for p in projects:
            if p.kind in ('biological', 'mixture') and (p.path / 'benchmark' / 'depth' / 'depth.tsv').is_file():
                f = depth_figure(p)
                if f:
                    out.append(f)
    return out
