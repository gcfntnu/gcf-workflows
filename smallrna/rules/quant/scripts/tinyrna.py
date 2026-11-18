#!/usr/bin/env python
"""tiny-count steps of the smallrna tinyrna route (#373).

features: write the Features Sheet (tiny-count's counting rules) from the route settings. One rule per 5' offset,
    all on mature miRNA features (GFF type miRNA, any attribute) with one hierarchy and classifier:
    `5' anchored, <offset>, <extension>` moves the feature's 5' end by the offset (positive: towards the 3' end) and
    its 3' end out by the extension, and selects alignments that start at the moved 5' end and end within the moved
    3' end. So with offsets -1,0,1 and extension 4, a read counts for the miRNA when it starts within one base of the
    annotated 5' end and ends at most 4 bases beyond the annotated 3' end (reads ending earlier count too). Strand
    and length are the other selectors; mismatches are left to the alignment.
count: run tiny-count standalone on one sample's alignments and keep its tables under fixed names. tiny-count writes
    into tiny-count_<timestamp>/ with the timestamp in every file name, needs every Paths File key (also those only
    the full tinyRNA pipeline uses) and exits 0 after an error, so its log and outputs are checked here.
join: join the samples' feature_counts.csv into one table (Feature ID, Classifier, Feature Name, one column per
    sample).

Usage: tinyrna.py features --offsets -1,0,1 --extension 4 --length 16-28 --strand sense > features.csv
       tinyrna.py count --sample S --alignments BAM --gff GFF3 --features CSV [--normalize-by-genomic-hits T/F]
                        [--normalize-by-feature-hits T/F] --output DIR
       tinyrna.py join --samples S1 S2 --tables T1 T2 > tinyrna_counts.tsv
"""
import argparse
import csv
from pathlib import Path
import re
import shutil
import subprocess
import sys

FEATURES_HEADER = ('Select for...', 'with value...', 'Classify as...', 'Source Filter', 'Type Filter', 'Hierarchy',
                   'Overlap', 'Mismatches', 'Strand', "5' End Nucleotide", 'Length')
STRANDS = ('sense', 'antisense', 'both')
# tiny-count's run tables (<prefix>_<timestamp>_<name>.csv) and per-library tables (..._<library>_<name>.csv)
RUN_TABLES = ('feature_counts', 'counts_by_rule', 'summary_stats', 'alignment_stats')
LIBRARY_TABLES = ('assigned_nt_len_dist', 'mapped_nt_len_dist')
PREFIX = 'tinyrna'


def fail(message):
    sys.exit('tinyrna.py: ' + message)


def parse_offsets(text):
    try:
        offsets = [int(v) for v in str(text).split(',') if v.strip()]
    except ValueError:
        fail("5' offsets must be integers separated by commas, got {!r}".format(text))
    if not offsets or len(set(offsets)) != len(offsets):
        fail("5' offsets must be one or more distinct integers, got {!r}".format(text))
    return sorted(offsets)


def features_sheet(offsets, extension, length, strand):
    """Rows of the Features Sheet (header first)."""
    if int(extension) < 0:
        fail("the 3' extension must be 0 or more, got {}".format(extension))
    if strand not in STRANDS:
        fail('strand must be one of {}, got {!r}'.format(', '.join(STRANDS), strand))
    if not re.fullmatch(r'\s*(all|\d+(\s*-\s*\d+)?(\s*,\s*\d+(\s*-\s*\d+)?)*)\s*', str(length)):
        fail('length must be all, a length, a range or a list (e.g. 16-28), got {!r}'.format(length))
    rows = [FEATURES_HEADER]
    for offset in parse_offsets(offsets):
        rows.append(('all', 'all', 'miRNA', '', 'miRNA', '1', "5' anchored, {}, {}".format(offset, int(extension)), '',
                     strand, 'all', str(length).replace(' ', '')))
    return rows


def paths_file(samples_csv, features_csv, gff):
    """tiny-count's Paths File; the keys of the full pipeline are present but empty, as tiny-count requires them."""
    return ('samples_csv: {}\nfeatures_csv: {}\ngff_files:\n- path: {}\n  alias: [Name]\n'
            'run_directory: run_directory\ntmp_directory:\nebwt:\nreference_genome_files:\nadapter_fasta:\n'
            'plot_style_sheet:\n').format(samples_csv, features_csv, gff)


def library_name(sample):
    """tiny-count's library name: Sample/Group Name and Replicate Number joined by _rep_."""
    return '{}_rep_1'.format(sample)


def collect(run_dir, sample, output):
    """Move tiny-count's tables from run_dir (tiny-count_<timestamp>/) to output under fixed names, and the config
    copies tiny-count made to output/config (all are rule outputs). Returns the moved table names."""
    found = {}
    for path in run_dir.glob('{}_*.csv'.format(PREFIX)):
        for name in RUN_TABLES + LIBRARY_TABLES:
            if path.name.endswith('_{}.csv'.format(name)):
                if name in LIBRARY_TABLES and not path.name.endswith('_{}_{}.csv'.format(library_name(sample), name)):
                    continue
                found[name] = path
    missing = [n for n in RUN_TABLES + LIBRARY_TABLES if n not in found]
    if missing:
        fail('tiny-count wrote no {} table(s) in {}'.format(', '.join(missing), run_dir))
    if not (run_dir / 'config').is_dir():
        fail('tiny-count wrote no config folder in {}'.format(run_dir))
    for name, path in found.items():
        shutil.move(str(path), str(output / '{}.csv'.format(name)))
    if (output / 'config').exists():
        shutil.rmtree(output / 'config')
    shutil.move(str(run_dir / 'config'), str(output / 'config'))
    return sorted(found)


def count(args):
    output = args.output.resolve()
    work = output / 'work'
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    alignments = args.alignments.resolve()
    with open(work / 'samples.csv', 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(('Input Files', 'Sample/Group Name', 'Replicate Number', 'Control', 'Normalization'))
        writer.writerow((str(alignments), args.sample, '1', '', ''))
    shutil.copy(args.features, work / 'features.csv')
    (work / 'paths.yml').write_text(paths_file('samples.csv', 'features.csv', args.gff.resolve()))
    log = output / 'tiny-count.log'
    command = ['tiny-count', '-pf', 'paths.yml', '-o', PREFIX, '-ng', args.normalize_by_genomic_hits,
               '-nf', args.normalize_by_feature_hits]
    with open(log, 'w') as fh:
        fh.write('# {}\n'.format(' '.join(command)))
        fh.flush()
        result = subprocess.run(command, cwd=work, stdout=fh, stderr=subprocess.STDOUT)
    text = log.read_text()
    if result.returncode != 0 or 'Traceback (most recent call last)' in text:
        fail('tiny-count failed for {} (exit {}); see {}'.format(args.sample, result.returncode, log))
    runs = sorted(work.glob('tiny-count_*'))
    if len(runs) != 1:
        fail('expected one tiny-count_<timestamp> directory in {}, found {}'.format(work, len(runs)))
    collect(runs[0], args.sample, output)
    shutil.rmtree(work)


def read_counts(path, sample):
    """{(Feature ID, Classifier): (Feature Name, count)} from one sample's feature_counts.csv."""
    with open(path, newline='') as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if header[:3] != ['Feature ID', 'Classifier', 'Feature Name'] or len(header) != 4:
            fail('{}: expected Feature ID, Classifier, Feature Name and one library column, got {}'.format(path, header))
        if header[3] != library_name(sample):
            fail('{}: library column {} is not sample {}'.format(path, header[3], sample))
        return {(row[0], row[1]): (row[2], row[3]) for row in reader if row}


def join(args):
    if len(args.samples) != len(args.tables):
        fail('give one table per sample')
    tables = [read_counts(t, s) for s, t in zip(args.samples, args.tables)]
    keys = sorted(set().union(*tables))
    writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    writer.writerow(['Feature ID', 'Classifier', 'Feature Name'] + list(args.samples))
    for key in keys:
        name = next((t[key][0] for t in tables if key in t), '')
        writer.writerow(list(key) + [name] + [t[key][1] if key in t else '0' for t in tables])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('features')
    p.add_argument('--offsets', required=True)
    p.add_argument('--extension', required=True, type=int)
    p.add_argument('--length', required=True)
    p.add_argument('--strand', required=True)
    p = sub.add_parser('count')
    p.add_argument('--sample', required=True)
    p.add_argument('--alignments', required=True, type=Path)
    p.add_argument('--gff', required=True, type=Path)
    p.add_argument('--features', required=True, type=Path)
    p.add_argument('--normalize-by-genomic-hits', default='T', choices=('T', 'F'))
    p.add_argument('--normalize-by-feature-hits', default='T', choices=('T', 'F'))
    p.add_argument('--output', required=True, type=Path)
    p = sub.add_parser('join')
    p.add_argument('--samples', required=True, nargs='+')
    p.add_argument('--tables', required=True, nargs='+', type=Path)
    args = parser.parse_args()
    if args.command == 'features':
        csv.writer(sys.stdout, lineterminator='\n').writerows(features_sheet(args.offsets, args.extension, args.length, args.strand))
    elif args.command == 'count':
        count(args)
    else:
        join(args)


if __name__ == '__main__':
    main()
