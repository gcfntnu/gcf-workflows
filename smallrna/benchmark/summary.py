"""Summary across evaluated benchmark projects: one section per dataset kind, a row per project and the kind's key
metrics per method, read from each project's `benchmark/metrics.tsv` and its configuration (`benchmark/dataset.yaml`,
`config.yaml`; lanes from `benchmark/samples.tsv`), with the reads its run processed as recorded by the evaluation.
Written as Markdown or as one self-contained HTML file, by the output's extension; nothing is written into the
projects.

Datasets have different methods, and a metric may be missing from a project's table: both give empty cells.
"""
from collections import deque
from contextlib import contextmanager
import csv
import html
import json
import os
from pathlib import Path
import re
import secrets
import statistics
import urllib.parse

import evaluation

import yaml

# (kind, title, key metrics, what they show); a pool's kind is its `pool`
SECTIONS = [
    ('synthetic', 'Synthetic data',
     ['recovered', 'precision', 'exact', 'misassigned', 'unplaced', 'near_spill', 'edger_called', 'edger_false_calls', 'edger_fdr',
      'isomir_recovered', 'decoy_as_mirna'],
     'Counts are compared per sample and miRBase accession. `recovered` sums the smaller of observed and true counts '
     'for each pair, divided by the total true count. `exact` is the fraction of positive-truth pairs with matching counts. '
     '`misassigned` and `unplaced` count reads on accessions absent from that sample’s truth and on names without a miRBase accession. '
     'The genome-wide dataset (#412) adds `precision` (the smaller counts over the observed ones), `near_spill` (reads above the truth '
     'on near-identical partners); there `misassigned` and `unplaced` are shares of the observed reads. Differential expression '
     'between its two groups (#447, `evaluate --de`) uses the tests of the ratiometric pools, so the two can be compared: with '
     'edgeR\'s exact test on the workflow\'s library sizes, `edger_called` (miRNAs designed at least 2-fold apart, BH-adjusted p at '
     'most 0.05 in the designed direction), `edger_false_calls` (miRNAs designed not to change and called) and `edger_fdr` (the '
     'share of calls on unchanged miRNAs or in the wrong direction). edgeR\'s GLM test, limma-voom, the tests adjusted for the '
     'nuisance factor and the built-in check `de_*` are in each report. `isomir_recovered` (#454): the share of its true isomiR '
     'reads (no substitutions) a method reports as their own sequence on their miRNA; empty for methods without isomiR output. '
     '`decoy_as_mirna` (#412 PR 2): in probe samples of reads of other RNA classes (tRNA, rRNA, Y RNA, snoRNA, other ncRNA, '
     'mRNA, piRNA, precursor), the share a method counts as miRNAs (control subtracted); each report has it per class and '
     'source, the miRNAs they land on, and the recovery of miRBase matures inside other-class molecules.'),
    ('equimolar', 'Equimolar pools', ['detected', 'cv', 'within_2fold', 'spearman_vs_oracle', 'outside', 'unplaced', 'shifted_within_2fold',
                                      'shifted_isomir_recovered'],
     'Scored `exact` entries share no miRNA accession with any other exact or shifted entry. The metrics are averaged across samples: '
     '`detected` and `cv` across these entries, `within_2fold` and `spearman_vs_oracle` against '
     'the oracle (the filtered reads counted by exact match to the pool sequences), over the entries the oracle counts in '
     'that sample only, so a shallow sample with few such entries can reach 1.0; `within_2fold` is the share with '
     '(method + 1) / (oracle + 1) within 2-fold. `outside` and `unplaced` are '
     'unweighted means of per-sample shares. Each denominator is all miRNA reads the method reported for that sample. '
     'The numerators are reads assigned to accessions with no exact or shifted entry in the analyzed pool and reads on names without a miRBase accession, respectively. '
     '`shifted_within_2fold`: the same for isomiRs of known amount (shifted entries), also with one added read. '
     '`shifted_isomir_recovered` (#454): the share of the oracle\'s reads of the shifted sequences a method reports as that '
     'sequence on its miRNA; empty for methods without isomiR output.'),
    ('ratiometric', 'Ratiometric pools', ['median_abs_error', 'spearman_vs_design', 'called', 'false_calls', 'shifted_called',
                                          'edger_called', 'edger_false_calls', 'edger_fdr', 'edger_fc_shift_vs_oracle',
                                          'edger_fc_error_vs_oracle', 'shifted_isomir_recovered'],
     'Observed log2(A/B) against the design: `median_abs_error`, `spearman_vs_design`, `called` (designed at least 2-fold '
     'apart and observed at least 2-fold in the designed direction), `false_calls` (designed less than 2-fold apart, '
     'observed at least 2-fold); `shifted_called` applies the same direction requirement to shifted sequences. '
     'Differential expression (#414, `evaluate --de`), with edgeR\'s exact test on the workflow\'s library sizes: `edger_called` '
     '(designed at least 2-fold apart, BH-adjusted p at most 0.05 in the designed direction), `edger_false_calls` (designed '
     '1:1 and called) and `edger_fdr` (the share of calls on 1:1 sequences or in the wrong direction); against the oracle\'s '
     'fits on the same runs (#444), `edger_fc_shift_vs_oracle` (the median difference from its log2 fold changes, from the '
     'library sizes) and `edger_fc_error_vs_oracle` (the median absolute deviation from that shift, from the counts). edgeR\'s '
     'GLM test, limma-voom and the tests within each lab are in each report. `shifted_isomir_recovered` as for equimolar pools.'),
    ('biological', 'Biological samples', ['mirna_reads', 'unplaced', 'detected', 'detected_all_samples', 'top10_share', 'isomir_share',
                                          'edger_null_called', 'voom_null_called',
                                          'knockout_mirna_share', 'knockout_surviving', 'detected_unflagged', 'read_in_other_class'],
     'No known composition: `mirna_reads` is the rounded median per sample of reads assigned to a miRBase accession. '
     'The other metrics describe the share on names without a miRBase accession, '
     'miRNAs detected (at least 10 reads; median over samples) and in every sample, and the share of the 10 most abundant. '
     'Agreement between methods is in each report. `isomir_share` (#454): the share of a method\'s miRNA reads that are not the '
     'mature sequence; empty for methods without isomiR output. With null comparisons (#449, `evaluate --de`): '
     '`edger_null_called` and `voom_null_called`, the median share of tested miRNAs that edgeR\'s exact test and limma-voom call '
     'when one group\'s samples, split in two at random, are tested against each other (every call is false). A dataset with '
     'knockout cells (#452): `knockout_mirna_share`, the miRNA reads per million filtered reads in the knockout that should '
     'remove miRNAs (DROSHA) relative to its control cells, and `knockout_surviving`, the miRNAs that keep more than 25% there. '
     'With the other RNA classes declared (#415): `detected_unflagged`, `detected` without the miRBase entries whose mature lies '
     'inside a tRNA, rRNA, Y RNA, snoRNA or other ncRNA, and `read_in_other_class`, the share of a method\'s miRNA reads (isomiR '
     'output) whose read lies inside such a molecule.'),
    ('mixture', 'Tissue mixtures', ['mirna_reads', 'detected', 'linearity_error', 'replicate_error', 'within_1_5fold', 'monotonic',
                                    'edger_glm_called', 'edger_glm_wrong_direction', 'edger_glm_called_shared'],
     'Total RNA of two tissues of each subject mixed at known fractions (#451). No known composition: each mixture is predicted '
     'from the method\'s own counts of the subject\'s two pure samples, with a content ratio fitted per method and subject, so the scores '
     'show linear and consistent counting, not correct assignment. `linearity_error` is the median |log2(observed / predicted)| '
     'over miRNA-sample pairs with at least 10 predicted reads and `within_1_5fold` the share within 1.5-fold; '
     '`replicate_error`, the same error between technical replicates, is its floor. `monotonic`: tissue-specific miRNAs whose '
     'subject series never move against the trend. Differential expression between the mixtures\' fractions, paired by subject '
     '(`evaluate --de`), with edgeR\'s GLM test: `edger_glm_called` (miRNAs the method expects at least 2-fold apart, called in '
     'that direction), `edger_glm_wrong_direction` (the share of its calls against the expected direction) and '
     '`edger_glm_called_shared` (`called` over the miRNAs every method expects at least 2-fold apart).'),
]
FIXED = ['project', 'dataset', 'kit', 'samples', 'reads processed', 'methods']


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
        raise ValueError('{}: dataset kind {!r} (pool {!r}) is not one of synthetic, pool (equimolar, ratiometric), biological, mixture'.format(
            project, spec.get('kind'), spec.get('pool')))
    table, record, intervals, routes = evaluation.read_evaluation(bench)
    columns = table[0][1:]
    values = {(row[0], c): v for row in table[1:] for c, v in zip(columns, row[1:]) if v != ''}
    for ci in intervals:  # bootstrap intervals (#413): the value with its interval
        if (ci['metric'], ci['row']) in values and ci['low'] != '' and ci['high'] != '':
            values[(ci['metric'], ci['row'])] = '{} ({}–{})'.format(values[(ci['metric'], ci['row'])], ci['low'], ci['high'])
    config = yaml.safe_load((project / 'config.yaml').read_text())
    configured_samples = config['samples']
    samples = len(configured_samples)
    if kind == 'synthetic' and (spec.get('synthetic') or {}).get('mode') == 'genome':  # #412: mixed and probe samples
        mixed = sum(1 for row in read_tsv(bench / 'samples.tsv') if row['kind'] == 'mixed')
        samples = '{} mixed, {} probe'.format(mixed, samples - mixed)
    if kind == 'biological':
        rows = read_tsv(bench / 'samples.tsv')
        if len(rows) != samples or {row['sample'] for row in rows} != set(configured_samples):
            raise ValueError(f'{project}: samples.tsv must have exactly one row per configured sample')
        libraries = {row['library'] for row in rows if row['lane'] != 'merged'}
        if libraries:
            samples = '{} lanes of {} samples'.format(samples, len(libraries))
    # declared routes not evaluated (#448): their state fills their cells; none for evaluations before routes.tsv
    status = {row['route']: row['state'] for row in routes if row['state'] != 'evaluated'}
    # one draw of a synthetic set (#455): the declaration without its seed and route selection identifies the set
    synthetic = dict(spec.get('synthetic') or {})
    seed = record.get('synthetic_seed', synthetic.pop('seed', None))
    synthetic.pop('seed', None)
    declaration = json.dumps(dict({k: v for k, v in spec.items() if k not in ('methods', 'not_selected')}, synthetic=synthetic),
                             sort_keys=True, default=str) if spec.get('kind') == 'synthetic' else None
    report = os.path.relpath((bench / 'report.md').resolve(), out_dir)
    return kind, {'project': project.resolve().name, 'link': urllib.parse.quote(report), 'dataset': spec['name'], 'kit': spec['libprepkit'],
                  'samples': samples, 'reads processed': '' if record.get('reads_processed') is None else '{:,}'.format(record['reads_processed']),
                  'methods': ', '.join(c for c in columns if c not in ('oracle', 'truth')), 'columns': columns, 'values': values,
                  'status': status, 'seed': seed, 'declaration': declaration}


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
        methods = list(dict.fromkeys([c for row in rows for c in row['columns']] + [m for row in rows for m in row['status']]))
        for row in rows:
            row['values'].update({(k, m): state for k in keys for m, state in row['status'].items()})
        metrics = [(k, [m for m in methods if any((k, m) in row['values'] for row in rows)]) for k in keys]
        metrics = [(k, ms) for k, ms in metrics if ms]
        sections.append((title, note, metrics, rows, over_seeds(rows, metrics)))
    return sections


def number(value):
    """The value of a metric cell without its interval, or None (empty, a route state)."""
    try:
        return float(str(value).split(' ')[0])
    except ValueError:
        return None


def over_seeds(rows, metrics):
    """Projects of one synthetic declaration that differ only by the seed of the draw (#455), two seeds or more:
    [(dataset, seeds, [(metric, {method: (mean, low, high)}, order)])], `order` how many pairs of the methods with a
    value in every seed keep their order (higher, lower or equal) in every seed, '' with fewer than two such methods.
    Two projects with the same seed stop it: each seed is one draw."""
    groups = {}
    for row in rows:
        if row['declaration'] is not None and row['seed'] is not None:
            groups.setdefault(row['declaration'], []).append(row)
    out = []
    for members in groups.values():
        seeds = sorted(row['seed'] for row in members)
        repeated = sorted({seed for seed in seeds if seeds.count(seed) > 1})
        if repeated:
            raise ValueError('projects {} draw synthetic set {} with the same seed {}; a seed is one draw, so give each project its own'.format(
                ', '.join(row['project'] for row in members if row['seed'] in repeated), members[0]['dataset'], ', '.join(map(str, repeated))))
        if len(seeds) < 2:
            continue
        lines = []
        for k, ms in metrics:
            values = {m: [number(row['values'].get((k, m))) for row in members] for m in ms}
            complete = {m: v for m, v in values.items() if v and None not in v}
            spread = {m: (statistics.mean(v), min(v), max(v)) for m, v in complete.items()}
            sign = lambda x: (x > 0) - (x < 0)
            pairs = [(a, b) for i, a in enumerate(complete) for b in list(complete)[i + 1:]]
            kept = sum(len({sign(complete[a][i] - complete[b][i]) for i in range(len(members))}) == 1 for a, b in pairs)
            lines.append((k, spread, '{} of {} pairs'.format(kept, len(pairs)) if pairs else ''))
        out.append((members[0]['dataset'], seeds, lines))
    return out


def seed_cell(spread):
    return '' if spread is None else '{:.4g} ({:.4g}–{:.4g})'.format(*spread)


INTRO = ('Key metrics of each evaluated project, from its `metrics.tsv`; its `report.md` has them all, and '
         '`smallrna/benchmark/README.md` defines them. `reads processed` is the reads the evaluated run filtered, from '
         'its cutadapt reports as recorded by `evaluate` (empty for another filter route, without the reports, or for an '
         'evaluation older than this record). `not run` or `stopped` in a metric cell: a route the dataset declares that '
         'was not selected or not run (`benchmark/not_run.tsv`), or whose checks failed; the report\'s Routes section has '
         'the reason (#448). An empty metric cell: the metric is not in the project\'s table, the method is not declared '
         'there, or the evaluation predates the route states. `oracle` is the pool reads counted by exact '
         'match, the reference for the methods; `truth` (genome-wide synthetic data) is the true counts analysed as a '
         'method, the best its design allows. A value with a range in parentheses has its 95% bootstrap interval, '
         'from the evaluation (`metrics_ci.tsv`; differences between methods are in each report).')


SEEDS_NOTE = ('The projects of this synthetic declaration that differ only by the seed of the draw (#455): per key metric '
              'and route, the mean over the seeds and their range in parentheses; `order`, how many pairs of the routes '
              'with a value in every seed keep their order (higher, lower or equal) in every seed.')


def markdown_cell(value):
    # Character references stay inside a GFM cell, including inside a link label.
    text = html.escape(str(value)).replace('\r', ' ').replace('\n', ' ')
    for character in ('\\', '|', '[', ']', '*', '`', '~'):
        text = text.replace(character, f'&#{ord(character)};')
    # Intraword underscores are literal in GFM. Escape delimiter runs at word boundaries.
    return re.sub(r'(?<![a-zA-Z0-9])_|_(?![a-zA-Z0-9])', '&#95;', text)


def markdown(sections):
    lines = ['# Small-RNA benchmark summary', '', INTRO]
    for title, note, metrics, rows, seeds in sections:
        header = FIXED + ['{}: {}'.format(k, m) for k, ms in metrics for m in ms]
        lines += ['', '## ' + title, '', note, '', '| ' + ' | '.join(map(markdown_cell, header)) + ' |', '|---' * len(header) + '|']
        for row in rows:
            cells = ['[{}]({})'.format(markdown_cell(row['project']), row['link'])] + [markdown_cell(row[k]) for k in FIXED[1:]]
            lines.append('| ' + ' | '.join(cells + [markdown_cell(row['values'].get((k, m), '')) for k, ms in metrics for m in ms]) + ' |')
        for dataset, used, table in seeds:
            methods = list(dict.fromkeys(m for _, spread, _ in table for m in spread))
            lines += ['', '### {} over seeds {}'.format(markdown_cell(dataset), ', '.join(map(str, used))), '', SEEDS_NOTE, '',
                      '| metric | ' + ' | '.join(map(markdown_cell, methods)) + ' | order |', '|---' * (len(methods) + 2) + '|']
            lines += ['| ' + ' | '.join([markdown_cell(k)] + [markdown_cell(seed_cell(spread.get(m))) for m in methods] + [order]) + ' |'
                      for k, spread, order in table]
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
    for title, note, metrics, rows, seeds in sections:
        top = ''.join('<th rowspan="2">{}</th>'.format(k) for k in FIXED)
        top += ''.join('<th colspan="{}">{}</th>'.format(len(ms), html.escape(k)) for k, ms in metrics)
        sub = ''.join('<th>{}</th>'.format(html.escape(m)) for _, ms in metrics for m in ms)
        trs = []
        for row in rows:
            cells = ['<td><a href="{}">{}</a></td>'.format(html.escape(row['link']), html.escape(row['project']))]
            cells += ['<td{}>{}</td>'.format(' class="num"' if k in ('samples', 'reads processed') else '', html.escape(str(row[k]))) for k in FIXED[1:]]
            cells += ['<td class="num">{}</td>'.format(html.escape(row['values'].get((k, m), ''))) for k, ms in metrics for m in ms]
            trs.append('<tr>{}</tr>'.format(''.join(cells)))
        body += ['<h2>{}</h2>'.format(html.escape(title)), '<p>{}</p>'.format(inline_code(note)),
                 '<div class="scroll"><table><thead><tr>{}</tr><tr>{}</tr></thead><tbody>'.format(top, sub)] + trs + ['</tbody></table></div>']
        for dataset, used, table in seeds:
            methods = list(dict.fromkeys(m for _, spread, _ in table for m in spread))
            head = ''.join('<th>{}</th>'.format(html.escape(m)) for m in ['metric'] + methods + ['order'])
            trs = ['<tr><td>{}</td>{}<td>{}</td></tr>'.format(html.escape(k), ''.join('<td class="num">{}</td>'.format(
                html.escape(seed_cell(spread.get(m)))) for m in methods), html.escape(order)) for k, spread, order in table]
            body += ['<h3>{} over seeds {}</h3>'.format(html.escape(dataset), html.escape(', '.join(map(str, used)))),
                     '<p>{}</p>'.format(inline_code(SEEDS_NOTE)),
                     '<div class="scroll"><table><thead><tr>{}</tr></thead><tbody>'.format(head)] + trs + ['</tbody></table></div>']
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
    projects, output = [Path(p).resolve(strict=True) for p in projects], Path(output)
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
