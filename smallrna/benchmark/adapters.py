"""Turn each quantification method's native output into one common table.

Common table: method, sample, mature_name (as the method reports it), mature_id (miRBase accession, empty when the
name cannot be placed) and count, one row per method, sample and reported miRNA. Counts are each method's total for
the mature miRNA, canonical reads and isomiRs together. Matching between methods and the truth is by accession, never
by name, because the methods name miRNAs differently:

- UNITAS drops the species prefix (`miR-192-5p`);
- miRge3 reports merge groups of identical sequences as `name1/name2`; such a count is given to the first member and
  the group is kept in `mature_name`;
- mirflowz names a miRNA after its locus (`mmu-miR-194-2-5p`) and reports isomiRs as separate rows
  (`name|5' shift|3' shift|CIGAR|length|sequence`); its own extended annotation maps names to accessions;
- TEsmall reports mature miRNAs by miRBase name (`ftype` miRNA); reads assigned to the precursor (`hairpin`, e.g.
  3' isomiRs that extend past the mature end) are not mature counts and are left out.

Adding a method means adding one function to ADAPTERS that reads its output directory.
"""
import csv
from pathlib import Path
import re

COMMON_COLUMNS = ('method', 'sample', 'mature_name', 'mature_id', 'count')


def read_matrix(path, sep):
    """Feature-by-sample count table -> {feature: {sample: count}}; empty cells and NA (mirflowz writes NA for a
    feature without reads in a sample) are 0."""
    with open(path, newline='') as fh:
        rows = list(csv.reader(fh, delimiter=sep))
    samples = rows[0][1:]
    return {row[0].strip(): {s: 0.0 if v in ('', 'NA') else float(v) for s, v in zip(samples, row[1:])} for row in rows[1:] if row}


def unitas(quant_dir, name_to_id, prefix, **_):
    rows = []
    for feature, counts in read_matrix(Path(quant_dir) / 'unitas' / 'mir_counts.tsv', '\t').items():
        name = feature if feature.startswith(prefix) else prefix + feature
        rows += [('unitas', s, feature, name_to_id.get(name, ''), c) for s, c in counts.items()]
    return rows


def mirge3(quant_dir, name_to_id, **_):
    """The workflow's miR.Counts.csv (built from miRge3's mapped.csv, #349), not miRge3's own miR.Counts.mirge3.csv."""
    rows = []
    for feature, counts in read_matrix(Path(quant_dir) / 'mirge3' / 'miR.Counts.csv', ',').items():
        if feature.lower() in ('mirna', 'mirna_name'):
            continue
        first = feature.split('/')[0]
        rows += [('mirge3', s, feature, name_to_id.get(first, ''), c) for s, c in counts.items()]
    return rows


def mirflowz_names(gff):
    """Name -> accession from mirflowz's extended miRNA annotation (Name includes the locus, Alias is MIMAT)."""
    mapping = {}
    with open(gff) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if not line.startswith('#') and len(fields) == 9:
                attributes = dict(a.split('=', 1) for a in fields[8].split(';') if '=' in a)
                if 'Name' in attributes and 'Alias' in attributes:
                    mapping[attributes['Name']] = attributes['Alias']
    return mapping


def mirflowz(quant_dir, name_to_id, mirflowz_gff=None, **_):
    locus_names = mirflowz_names(mirflowz_gff) if mirflowz_gff else {}
    totals = {}
    for feature, counts in read_matrix(Path(quant_dir) / 'mirflowz' / 'all_mirna_counts.tab', '\t').items():
        name = feature.split('|')[0]
        for sample, count in counts.items():
            totals[(name, sample)] = totals.get((name, sample), 0) + count
    return [('mirflowz', s, name, locus_names.get(name) or name_to_id.get(name, ''), c) for (name, s), c in totals.items()]


def tesmall(quant_dir, name_to_id, **_):
    """TEsmall's joined table (ftype, fid, samples); mature miRNAs only, by miRBase name."""
    rows = []
    with open(Path(quant_dir) / 'tesmall' / 'tesmall_counts.tsv', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        samples = next(reader)[2:]
        for row in reader:
            if row and row[0] == 'miRNA':
                rows += [('tesmall', s, row[1], name_to_id.get(row[1], ''), float(c or 0)) for s, c in zip(samples, row[2:])]
    return rows


ADAPTERS = {'unitas': unitas, 'mirge3': mirge3, 'mirflowz': mirflowz, 'tesmall': tesmall}


def common_table(methods, quant_dir, name_to_id, prefix, mirflowz_gff=None):
    unknown = [m for m in methods if m not in ADAPTERS]
    if unknown:
        raise ValueError('no benchmark adapter for method(s) {}; available: {}'.format(', '.join(unknown), ', '.join(ADAPTERS)))
    rows = []
    for method in methods:
        rows += ADAPTERS[method](quant_dir, name_to_id, prefix=prefix, mirflowz_gff=mirflowz_gff)
    return [row for row in rows if row[4]]
