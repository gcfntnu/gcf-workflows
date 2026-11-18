"""Summary across evaluated benchmark projects: one section per dataset kind, a row per project and the kind's key
metrics per method, read from each project's `benchmark/metrics.tsv` and its configuration (`benchmark/dataset.yaml`,
`config.yaml`; input reads from `benchmark/samples.tsv`, or `truth.tsv` for synthetic data). Written as Markdown or
as one self-contained HTML file, by the output's extension; nothing is written into the projects.

Datasets have different methods, and a metric may be missing from a project's table: both give empty cells.
"""
import csv
import html
import os
from pathlib import Path
import urllib.parse

import yaml

# (kind, title, key metrics, what they show); a pool's kind is its `pool`
SECTIONS = [
    ('synthetic', 'Synthetic data', ['recovered', 'exact', 'misassigned', 'unplaced'],
     'Against the truth per read: `recovered` (share of the true reads counted on the right miRNA), `exact` (share of '
     'miRNAs counted exactly), `misassigned` and `unplaced` (reads on miRNAs not in the truth, on names without a miRBase '
     'accession).'),
    ('equimolar', 'Equimolar pools', ['detected', 'cv', 'within_2fold', 'spearman_vs_oracle', 'outside', 'unplaced', 'shifted_within_2fold'],
     'Over the pool\'s `exact` sequences: `detected` and `cv` across the pool, `within_2fold` and `spearman_vs_oracle` against '
     'the oracle (the filtered reads counted by exact match to the pool sequences); `outside` and `unplaced`: shares of a '
     'method\'s miRNA reads on miRNAs of no pool sequence and on names without a miRBase accession; '
     '`shifted_within_2fold`: isomiRs of known amount against the oracle.'),
    ('ratiometric', 'Ratiometric pools', ['median_abs_error', 'spearman_vs_design', 'called', 'false_calls', 'shifted_called'],
     'Observed log2(A/B) against the design: `median_abs_error`, `spearman_vs_design`, `called` (designed at least 2-fold '
     'apart and observed so), `false_calls` (designed less than 2-fold apart, observed so); `shifted_called` for isomiRs.'),
    ('biological', 'Biological samples', ['mirna_reads', 'unplaced', 'detected', 'detected_all_samples', 'top10_share'],
     'No known composition, so agreement is not accuracy: reads on miRNAs, the share on names without a miRBase accession, '
     'miRNAs detected (at least 10 reads; median over samples) and in every sample, and the share of the 10 most abundant. '
     'Agreement between methods is in each report.'),
]
FIXED = ['project', 'dataset', 'kit', 'samples', 'reads', 'methods']


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def read_project(project, out_dir):
    """A project's row: its description (FIXED, with a link to its report relative to `out_dir`) and its metrics
    {(metric, method): value}, the empty ones left out."""
    bench = project / 'benchmark'
    if not (bench / 'dataset.yaml').is_file():
        raise ValueError('{}: no benchmark/dataset.yaml; not a benchmark project (`benchmark.py prepare`)'.format(project))
    if not (bench / 'metrics.tsv').is_file():
        raise ValueError('{}: no benchmark/metrics.tsv; run `benchmark.py evaluate --project {}`'.format(project, project))
    spec = yaml.safe_load((bench / 'dataset.yaml').read_text())
    kind = spec.get('pool') if spec.get('kind') == 'pool' else spec.get('kind')
    if kind not in [s[0] for s in SECTIONS]:
        raise ValueError('{}: dataset kind {!r} (pool {!r}) is not one of synthetic, pool (equimolar, ratiometric), biological'.format(
            project, spec.get('kind'), spec.get('pool')))
    with open(bench / 'metrics.tsv', newline='') as fh:
        table = list(csv.reader(fh, delimiter='\t'))
    columns = table[0][1:]
    values = {(row[0], c): v for row in table[1:] for c, v in zip(columns, row[1:]) if v != ''}
    samples = len(yaml.safe_load((project / 'config.yaml').read_text())['samples'])
    if kind == 'synthetic':
        reads = [row['count'] for row in read_tsv(bench / 'truth.tsv')]
    else:
        rows = read_tsv(bench / 'samples.tsv')
        reads = [row['reads'] for row in rows]
        libraries = {row['library'] for row in rows if row['lane'] != 'merged'} if kind == 'biological' else set()
        if libraries:
            samples = '{} lanes of {} samples'.format(samples, len(libraries))
    report = os.path.relpath(bench / 'report.md', out_dir)
    return kind, {'project': project.resolve().name, 'link': urllib.parse.quote(report), 'dataset': spec['name'], 'kit': spec['libprepkit'],
                  'samples': samples, 'reads': '{:,}'.format(sum(int(r) for r in reads)) if all(reads) else '',
                  'methods': ', '.join(c for c in columns if c != 'oracle'), 'columns': columns, 'values': values}


def summarize(projects, output):
    """[(title, note, [(metric, [method])], [row])]: the sections with projects, a metric's methods being those with a
    value in any of its rows (the oracle included), in order of appearance."""
    out_dir = output.resolve().parent
    names = [p.resolve().name for p in projects]
    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:
        raise ValueError('projects named alike: {}; the summary names a project by its directory'.format(', '.join(repeated)))
    inside = [str(p) for p in projects if output.resolve().is_relative_to(p.resolve())]
    if inside:
        raise ValueError('the output {} is inside the project {}; the summary writes nothing into the projects'.format(output, inside[0]))
    by_kind = {}
    for project in projects:
        kind, row = read_project(project, out_dir)
        by_kind.setdefault(kind, []).append(row)
    sections = []
    for kind, title, keys, note in SECTIONS:
        rows = by_kind.get(kind)
        if not rows:
            continue
        methods = list(dict.fromkeys(c for row in rows for c in row['columns']))
        metrics = [(k, [m for m in methods if any((k, m) in row['values'] for row in rows)]) for k in keys]
        sections.append((title, note, [(k, ms) for k, ms in metrics if ms], rows))
    return sections


INTRO = ('Key metrics of each evaluated project, from its `metrics.tsv`; its `report.md` has them all, and '
         '`smallrna/benchmark/README.md` defines them. An empty cell: the method was not run on that dataset, or the metric '
         'is not in its table. `oracle` is the pool reads counted by exact match, the reference for the methods.')


def markdown(sections):
    lines = ['# Small-RNA benchmark summary', '', INTRO]
    for title, note, metrics, rows in sections:
        header = FIXED + ['{}: {}'.format(k, m) for k, ms in metrics for m in ms]
        lines += ['', '## ' + title, '', note, '', '| ' + ' | '.join(header) + ' |', '|---' * len(header) + '|']
        for row in rows:
            cells = ['[{}]({})'.format(row['project'], row['link'])] + [str(row[k]) for k in FIXED[1:]]
            lines.append('| ' + ' | '.join(cells + [row['values'].get((k, m), '') for k, ms in metrics for m in ms]) + ' |')
    return '\n'.join(lines) + '\n'


STYLE = '''
:root { --bg: #ffffff; --fg: #1f2328; --muted: #59636e; --line: #d1d9e0; --head: #f6f8fa; --link: #0969da; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #0d1117; --fg: #e6edf3; --muted: #9198a1; --line: #3d444d; --head: #151b23; --link: #4493f8; }
}
body { background: var(--bg); color: var(--fg); font: 14px/1.5 system-ui, sans-serif; margin: 0 auto; padding: 16px; max-width: 1600px; }
h1 { font-size: 1.5em; } h2 { font-size: 1.2em; margin-top: 2em; }
p { color: var(--muted); max-width: 70em; }
a { color: var(--link); }
code { font-size: 0.95em; }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
th, td { border: 1px solid var(--line); padding: 3px 8px; white-space: nowrap; }
thead th { background: var(--head); }
td.num { text-align: right; }
'''


def inline_code(text):
    """Escaped text with `code` spans as <code>."""
    parts = html.escape(text).split('`')
    return ''.join('<code>{}</code>'.format(p) if i % 2 else p for i, p in enumerate(parts))


def html_page(sections):
    body = ['<h1>Small-RNA benchmark summary</h1>', '<p>{}</p>'.format(inline_code(INTRO))]
    for title, note, metrics, rows in sections:
        top = ''.join('<th rowspan="2">{}</th>'.format(k) for k in FIXED)
        top += ''.join('<th colspan="{}">{}</th>'.format(len(ms), html.escape(k)) for k, ms in metrics)
        sub = ''.join('<th>{}</th>'.format(html.escape(m)) for _, ms in metrics for m in ms)
        trs = []
        for row in rows:
            cells = ['<td><a href="{}">{}</a></td>'.format(html.escape(row['link']), html.escape(row['project']))]
            cells += ['<td{}>{}</td>'.format(' class="num"' if k in ('samples', 'reads') else '', html.escape(str(row[k]))) for k in FIXED[1:]]
            cells += ['<td class="num">{}</td>'.format(html.escape(row['values'].get((k, m), ''))) for k, ms in metrics for m in ms]
            trs.append('<tr>{}</tr>'.format(''.join(cells)))
        body += ['<h2>{}</h2>'.format(html.escape(title)), '<p>{}</p>'.format(inline_code(note)),
                 '<div class="scroll"><table><thead><tr>{}</tr><tr>{}</tr></thead><tbody>'.format(top, sub)] + trs + ['</tbody></table></div>']
    return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n<title>Small-RNA benchmark summary</title>\n'
            '<style>{}</style>\n</head>\n<body>\n{}\n</body>\n</html>\n'.format(STYLE, '\n'.join(body)))


def write(projects, output):
    """Write the summary of `projects` to `output` (.md or .html)."""
    projects, output = [Path(p) for p in projects], Path(output)
    render = {'.md': markdown, '.html': html_page}.get(output.suffix)
    if render is None:
        raise ValueError('--output {}: the extension must be .md or .html'.format(output))
    sections = summarize(projects, output)
    output.write_text(render(sections))
