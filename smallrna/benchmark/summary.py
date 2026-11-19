"""Summary across evaluated benchmark projects: one section per dataset kind, a row per project and the kind's key
metrics per method, read from each project's `benchmark/metrics.tsv` and its configuration (`benchmark/dataset.yaml`,
`config.yaml`; input reads from `benchmark/samples.tsv`, or `truth.tsv` for synthetic data). Written as Markdown or
as one self-contained HTML file, by the output's extension; nothing is written into the projects.

Datasets have different methods, and a metric may be missing from a project's table: both give empty cells.
"""
from collections import deque
from contextlib import contextmanager
import csv
import html
import os
from pathlib import Path
import re
import secrets
import urllib.parse

import evaluation

import yaml

# (kind, title, key metrics, what they show); a pool's kind is its `pool`
SECTIONS = [
    ('synthetic', 'Synthetic data', ['recovered', 'exact', 'misassigned', 'unplaced'],
     'Counts are compared per sample and miRBase accession. `recovered` sums the smaller of observed and true counts '
     'for each pair, divided by the total true count. `exact` is the fraction of positive-truth pairs with matching counts. '
     '`misassigned` and `unplaced` count reads on accessions absent from that sample’s truth and on names without a miRBase accession.'),
    ('equimolar', 'Equimolar pools', ['detected', 'cv', 'within_2fold', 'spearman_vs_oracle', 'outside', 'unplaced', 'shifted_within_2fold'],
     'Over the pool\'s `exact` sequences: `detected` and `cv` across the pool, `within_2fold` and `spearman_vs_oracle` against '
     'the oracle (the filtered reads counted by exact match to the pool sequences). `outside` and `unplaced` are '
     'unweighted means of per-sample shares. Each denominator is all miRNA reads the method reported for that sample. '
     'The numerators are reads assigned to accessions with no exact or shifted entry in the analyzed pool and reads on names without a miRBase accession, respectively. '
     '`shifted_within_2fold`: isomiRs of known amount against the oracle.'),
    ('ratiometric', 'Ratiometric pools', ['median_abs_error', 'spearman_vs_design', 'called', 'false_calls', 'shifted_called'],
     'Observed log2(A/B) against the design: `median_abs_error`, `spearman_vs_design`, `called` (designed at least 2-fold '
     'apart and observed at least 2-fold in the designed direction), `false_calls` (designed less than 2-fold apart, '
     'observed at least 2-fold); `shifted_called` applies the same direction requirement to shifted sequences.'),
    ('biological', 'Biological samples', ['mirna_reads', 'unplaced', 'detected', 'detected_all_samples', 'top10_share'],
     'No known composition: `mirna_reads` is the rounded median per sample of reads assigned to a miRBase accession. '
     'The other metrics describe the share on names without a miRBase accession, '
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
    spec = yaml.safe_load((bench / 'dataset.yaml').read_text())
    kind = spec.get('pool') if spec.get('kind') == 'pool' else spec.get('kind')
    if kind not in [s[0] for s in SECTIONS]:
        raise ValueError('{}: dataset kind {!r} (pool {!r}) is not one of synthetic, pool (equimolar, ratiometric), biological'.format(
            project, spec.get('kind'), spec.get('pool')))
    table = evaluation.read_metrics(bench)
    columns = table[0][1:]
    values = {(row[0], c): v for row in table[1:] for c, v in zip(columns, row[1:]) if v != ''}
    configured_samples = yaml.safe_load((project / 'config.yaml').read_text())['samples']
    samples = len(configured_samples)
    if kind == 'synthetic':
        rows = read_tsv(bench / 'truth.tsv')
        if {row['sample'] for row in rows} != set(configured_samples):
            raise ValueError(f'{project}: truth.tsv must cover every configured sample and no others')
        reads = [row['count'] for row in rows]
    else:
        rows = read_tsv(bench / 'samples.tsv')
        if len(rows) != samples or {row['sample'] for row in rows} != set(configured_samples):
            raise ValueError(f'{project}: samples.tsv must have exactly one row per configured sample')
        reads = [row['reads'] if row.get('reads_source') in ('copied', 'verified_declaration') else '' for row in rows]
        libraries = {row['library'] for row in rows if row['lane'] != 'merged'} if kind == 'biological' else set()
        if libraries:
            samples = '{} lanes of {} samples'.format(samples, len(libraries))
    report = os.path.relpath((bench / 'report.md').resolve(), out_dir)
    return kind, {'project': project.resolve().name, 'link': urllib.parse.quote(report), 'dataset': spec['name'], 'kit': spec['libprepkit'],
                  'samples': samples, 'reads': '{:,}'.format(sum(int(r) for r in reads)) if reads and all(reads) else '',
                  'methods': ', '.join(c for c in columns if c != 'oracle'), 'columns': columns, 'values': values}


def resolve_with_links(path):
    """Resolve a path and retain the traversed link entries, including links in link targets."""
    path = path.absolute()
    current = Path(path.anchor).resolve()
    remaining = deque(path.parts[1:])
    links = []
    while remaining:
        component = remaining.popleft()
        if component == '..':
            current = current.parent
            continue
        current /= component
        if current.is_symlink():
            links.append(current)
            if len(links) > 40:  # Linux pathname resolution limit
                raise ValueError(f'too many symbolic links in {path}')
            target = current.readlink()
            current = Path(target.anchor).resolve() if target.is_absolute() else current.parent
            remaining.extendleft(reversed(target.parts[1:] if target.is_absolute() else target.parts))
    return current, links


def summarize(projects, out_dir):
    """Sections with projects, metrics and method columns in order of appearance."""
    names = [p.resolve().name for p in projects]
    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:
        raise ValueError('projects named alike: {}; the summary names a project by its directory'.format(', '.join(repeated)))
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


def markdown_cell(value):
    # Character references stay inside a GFM cell, including inside a link label.
    text = html.escape(str(value)).replace('\r', ' ').replace('\n', ' ')
    for character in ('\\', '|', '[', ']', '*', '`', '~'):
        text = text.replace(character, f'&#{ord(character)};')
    # Intraword underscores are literal in GFM. Escape delimiter runs at word boundaries.
    return re.sub(r'(?<![a-zA-Z0-9])_|_(?![a-zA-Z0-9])', '&#95;', text)


def markdown(sections):
    lines = ['# Small-RNA benchmark summary', '', INTRO]
    for title, note, metrics, rows in sections:
        header = FIXED + ['{}: {}'.format(k, m) for k, ms in metrics for m in ms]
        lines += ['', '## ' + title, '', note, '', '| ' + ' | '.join(map(markdown_cell, header)) + ' |', '|---' * len(header) + '|']
        for row in rows:
            cells = ['[{}]({})'.format(markdown_cell(row['project']), row['link'])] + [markdown_cell(row[k]) for k in FIXED[1:]]
            lines.append('| ' + ' | '.join(cells + [markdown_cell(row['values'].get((k, m), '')) for k, ms in metrics for m in ms]) + ' |')
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


def directory_ancestors(directory):
    """Directory identities through the mounted tree, starting from the open output directory."""
    identities = set()
    current = os.dup(directory)
    try:
        while True:
            info = os.fstat(current)
            identity = (info.st_dev, info.st_ino)
            if identity in identities:  # At the filesystem namespace root, .. is the same directory.
                return identities
            identities.add(identity)
            parent = os.open('..', getattr(os, 'O_PATH', os.O_RDONLY) | os.O_DIRECTORY, dir_fd=current)
            os.close(current)
            current = parent
    finally:
        os.close(current)


def validate_destination(projects, entry, directory, directory_links):
    """Check the bound directory entry against the project locations and their link chains."""
    opened = os.fstat(directory)
    check_directory_identity(directory, entry.parent)
    ancestors = directory_ancestors(directory)
    final, file_links = resolve_with_links(entry)
    destinations = {entry, final, *directory_links, *file_links}
    for project in projects:
        bench = project / 'benchmark'
        project_dirs = {project.resolve(), bench.resolve()}
        # A mounted project may have a second absolute path. Check the opened
        # directory and its parents even for a previously unused output name.
        for root in project_dirs:
            try:
                info = root.stat()
            except FileNotFoundError:
                continue  # read_project diagnoses unprepared projects.
            if (info.st_dev, info.st_ino) in ancestors:
                raise ValueError(f'the output {entry} is inside the project {project}')
        sources = [project / 'config.yaml', *(bench / name for name in
                   ('dataset.yaml', 'metrics.tsv', 'report.md', '.evaluation.lock', 'evaluation.json', 'truth.tsv', 'samples.tsv'))]
        inputs = set()
        for source in sources:
            target, links = resolve_with_links(source)
            inputs.update([target, *links])
        # Compare directory entries, not file inodes. Replacing a separate hard
        # link outside a project preserves the project file.
        same_entry = any(path.name == entry.name and os.path.samestat(path.parent.stat(), opened) for path in inputs)
        if same_entry or destinations & inputs or any(path.is_relative_to(root) for path in destinations for root in project_dirs):
            raise ValueError(f'the output {entry} is inside the project {project} or aliases one of its inputs')


@contextmanager
def open_directory(path):
    """Pin a resolved Linux directory, refusing links substituted during traversal."""
    flags = getattr(os, 'O_PATH', os.O_RDONLY) | os.O_DIRECTORY | os.O_NOFOLLOW
    fd = os.open(path.anchor, flags)
    try:
        for part in path.parts[1:]:
            child = os.open(part, flags, dir_fd=fd)
            os.close(fd)
            fd = child
        yield fd
    finally:
        os.close(fd)


def check_directory_identity(directory, path, *, follow_symlinks=False):
    if not os.path.samestat(os.fstat(directory), path.stat(follow_symlinks=follow_symlinks)):
        raise ValueError(f'the summary directory moved during publication: {path}')


def write(projects, output):
    """Write the summary of `projects` to `output` (.md or .html)."""
    # Inspect the supplied filename before Path removes trailing slashes or ".".
    name = os.path.basename(output)
    render = {'.md': markdown, '.html': html_page}.get(Path(name).suffix)
    if render is None:
        raise ValueError('--output {}: the extension must be .md or .html'.format(output))
    projects, output = [Path(p) for p in projects], Path(output)
    requested_parent = output.parent
    parent, links = resolve_with_links(requested_parent)
    output = parent / output.name
    with open_directory(parent) as directory:
        # Kernel lookup of the original spelling must reach this directory.
        # Lexically removing ".." can otherwise skip a missing or non-directory entry.
        check_directory_identity(directory, requested_parent, follow_symlinks=True)
        validate_destination(projects, output, directory, links)
        sections = summarize(projects, output.parent)
        text = render(sections)
        check_directory_identity(directory, output.parent)
        temporary = '.summary-' + secrets.token_hex(16)
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666, dir_fd=directory)
        try:
            with os.fdopen(fd, 'w', encoding='utf-8') as stream:
                stream.write(text)
            validate_destination(projects, output, directory, links)
            os.replace(temporary, output.name, src_dir_fd=directory, dst_dir_fd=directory)
        except BaseException:
            os.unlink(temporary, dir_fd=directory)
            raise
        check_directory_identity(directory, output.parent)
    return output
