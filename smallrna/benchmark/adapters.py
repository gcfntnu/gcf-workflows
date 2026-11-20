"""Turn each quantification method's native output into one common table.

Common table: method, sample, mature_name (as the method reports it), mature_id (miRBase accession, empty when the
name cannot be placed) and count, one row per method, sample and reported miRNA. Counts are each method's total for
the mature miRNA, canonical reads and isomiRs together. Matching between methods and the truth is by accession, never
by name, because the methods name miRNAs differently:

- UNITAS drops the species prefix (`miR-192-5p`) and names reads it placed on a precursor rather than on a mature
  miRNA (with a modification, e.g. a substituted first base) after the hairpin and arm (`miR-16-1-5p`); such a name
  gets the organism's mature miRNA on that arm of the hairpin in the UNITAS reference (unitas_hairpins);
- miRge3 reports merge groups of identical sequences as `name1/name2`; such a count is given to the first member and
  the group is kept in `mature_name`;
- mirflowz names a miRNA after its locus (`mmu-miR-194-2-5p`) and reports isomiRs as separate rows
  (`name|5' shift|3' shift|CIGAR|length|sequence`); its own extended annotation maps names to accessions;
- TEsmall reports mature miRNAs by miRBase name (`ftype` miRNA); reads assigned to the precursor (`hairpin`, e.g.
  3' isomiRs that extend past the mature end) are not mature counts and are left out.
- tiny-count reports miRBase GFF feature IDs, one row per locus: the accession, with `_1`, `_2` for further loci of
  the same mature miRNA (`MIMAT0000416_1`); the loci are summed per accession;
- nf-core/smrnaseq's mirtop counting (`mirna.tsv`, smallrna/rules/quant/smrnaseq.smk) reports miRNAs by miRBase name,
  isomiRs summed.
- htseq (the genome counter of an earlier in-house pipeline, `htseq_counts.tsv`, smallrna/rules/quant/htseq.smk)
  reports miRNAs by miRBase name, isomiRs summed; a read counts in full for every name it overlaps.

Adding a method means adding one function to ADAPTERS that reads its output directory, and its count table to TABLES:
before scoring, benchmark.py checks that table's samples and the miRBase release in the method's reference.tsv (#407).
"""
import csv
import hashlib
import math
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


def unitas_hairpins(refdump):
    """The organism's hairpins in a UNITAS miRBase reference dump (`miRBase_hairpin_<Species>.fas`, names without the
    species prefix plus `_stem-loop`): {hairpin name: DNA sequence}."""
    found = [p for p in Path(refdump).glob('miRBase_hairpin_*.fas') if not p.name.startswith('miRBase_hairpin_non-')]
    if len(found) != 1:
        raise ValueError('{}: expected one miRBase_hairpin_<Species>.fas, found {}'.format(refdump, len(found)))
    hairpins, name = {}, None
    with open(found[0]) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0].replace('_stem-loop', '')
                hairpins[name] = ''
            elif name:
                hairpins[name] += line.strip().upper().replace('U', 'T')
    return hairpins


def unitas_place(feature, name_to_id, prefix, hairpins, mature):
    """(accession, step) for a UNITAS miRNA name: step `name` when it is a mature miRNA of the organism, `hairpin`
    when it is a hairpin and arm (`miR-16-1-5p`) and one of the organism's mature miRNAs ({name: (accession, DNA
    sequence)}) on that arm of the hairpin is the hairpin family's (hsa-miR-16-5p for miR-16-1) or the only one,
    otherwise ('', ''): no such hairpin, no mature miRNA on that arm, or several of other families."""
    name = feature if feature.startswith(prefix) else prefix + feature
    if name in name_to_id:
        return name_to_id[name], 'name'
    m = re.fullmatch(r'(?:{})?(.+)-([35]p)'.format(re.escape(prefix)), feature)
    hairpin = hairpins.get(m.group(1)) if m else None
    if not hairpin:
        return '', ''
    on_arm = {}
    for mature_name, (accession, seq) in mature.items():
        start = hairpin.find(seq)
        if start >= 0 and ('5p' if start + len(seq) / 2 < len(hairpin) / 2 else '3p') == m.group(2):
            on_arm[mature_name] = accession
    family = '{}{}-{}'.format(prefix, re.sub(r'-\d+$', '', m.group(1)), m.group(2))
    if family in on_arm:
        return on_arm[family], 'hairpin'
    return (next(iter(on_arm.values())), 'hairpin') if len(set(on_arm.values())) == 1 else ('', '')


def unitas(quant_dir, name_to_id, prefix, unitas_hairpins=None, mature=None, **_):
    """UNITAS's mir_counts.tsv; without hairpins (unitas_hairpins) and mature sequences, names are matched directly only."""
    rows = []
    for feature, counts in read_matrix(Path(quant_dir) / 'unitas' / 'mir_counts.tsv', '\t').items():
        accession, _ = unitas_place(feature, name_to_id, prefix, unitas_hairpins or {}, mature or {})
        rows += [('unitas', s, feature, accession, c) for s, c in counts.items()]
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


def tinyrna(quant_dir, name_to_id, **_):
    """The workflow's tinyrna_counts.tsv (tiny-count feature counts joined over samples); loci summed per accession."""
    shares, names = {}, {}
    with open(Path(quant_dir) / 'tinyrna' / 'tinyrna_counts.tsv', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        samples = next(reader)[3:]
        for row in reader:
            if not row or row[1] != 'miRNA':
                continue
            accession = re.sub(r'_\d+$', '', row[0])
            names.setdefault(accession, row[2].split(', ')[0])
            for sample, value in zip(samples, row[3:]):
                shares.setdefault((accession, sample), []).append(float(value or 0))
    return [('tinyrna', s, names[a], a, math.fsum(c)) for (a, s), c in shares.items()]


def smrnaseq(quant_dir, name_to_id, **_):
    """nf-core/smrnaseq's mirtop table as the route writes it (miRNA, samples), by miRBase name."""
    return named_counts('smrnaseq', Path(quant_dir) / 'smrnaseq' / 'mirna.tsv', name_to_id)


def htseq(quant_dir, name_to_id, **_):
    """The htseq route's joined mature counts (miRNA, samples), by miRBase name."""
    return named_counts('htseq', Path(quant_dir) / 'htseq' / 'htseq_counts.tsv', name_to_id)


def named_counts(method, path, name_to_id):
    """Common rows of a table of miRNA names and one column per sample. A header without unique sample columns, or a row
    whose field count differs from the header (e.g. a truncated file), stops the evaluation."""
    rows = []
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, None)
        samples = header[1:] if header and header[0] == 'miRNA' else []
        if not samples or len(set(samples)) != len(samples) or not all(samples):
            raise ValueError('{}: expected a miRNA column and unique sample columns, got {}'.format(path, header))
        for i, row in enumerate(reader, 2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError('{}:{}: {} fields, expected {}'.format(path, i, len(row), len(header)))
            rows += [(method, s, row[0], name_to_id.get(row[0], ''), float(c)) for s, c in zip(samples, row[1:])]
    return rows


ADAPTERS = {'unitas': unitas, 'mirge3': mirge3, 'mirflowz': mirflowz, 'tesmall': tesmall, 'tinyrna': tinyrna, 'smrnaseq': smrnaseq,
            'htseq': htseq}

# Each method's count table under the project's quant directory: (path, delimiter, columns before the samples).
TABLES = {'unitas': ('unitas/mir_counts.tsv', '\t', 1), 'mirge3': ('mirge3/miR.Counts.csv', ',', 1),
          'mirflowz': ('mirflowz/all_mirna_counts.tab', '\t', 1), 'tesmall': ('tesmall/tesmall_counts.tsv', '\t', 2),
          'tinyrna': ('tinyrna/tinyrna_counts.tsv', '\t', 3), 'smrnaseq': ('smrnaseq/mirna.tsv', '\t', 1),
          'htseq': ('htseq/htseq_counts.tsv', '\t', 1)}
# Methods whose reference.tsv is written by quant_reference after their count table, with that table's SHA-256.
TABLE_BOUND_RECORDS = ('mirflowz', 'tesmall', 'tinyrna', 'smrnaseq', 'htseq')


def table_samples(method, quant_dir):
    """The sample columns of a method's count table, in their order."""
    path, sep, skip = TABLES[method]
    with open(Path(quant_dir) / path, newline='') as fh:
        header = next(csv.reader(fh, delimiter=sep), [])
    return header[skip:]  # as written: the adapters score under the literal column names


def method_record(method, quant_dir, names):
    """{name: value} for the named fields of a method's quant/<method>/reference.tsv (written by unitas_reference,
    mirge3_reference and quant_reference): unique field and value columns, rows as wide as the header, and exactly one
    non-empty row per name; ValueError otherwise."""
    path = Path(quant_dir) / method / 'reference.tsv'
    with open(path, newline='') as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    header = rows[0] if rows else []
    if len(set(header)) != len(header) or 'field' not in header or 'value' not in header:
        raise ValueError('{}: expected unique field and value columns, got {}'.format(path, header))
    for i, row in enumerate(rows[1:], 2):
        if row and len(row) != len(header):
            raise ValueError('{}:{}: {} fields, expected {}'.format(path, i, len(row), len(header)))
    f, v = header.index('field'), header.index('value')
    out = {}
    for name in names:
        found = [row[v] for row in rows[1:] if row and row[f] == name]
        if len(found) != 1 or not found[0].strip():
            raise ValueError('{}: expected one non-empty {} row, found {}'.format(path, name, len(found)))
        out[name] = found[0]
    return out


def method_mirbase(method, quant_dir):
    """(miRNA database, release) a method annotated with, from its reference.tsv."""
    record = method_record(method, quant_dir, ('mirna_database', 'mirna_release'))
    return record['mirna_database'], record['mirna_release']


def table_sha256(method, quant_dir):
    """SHA-256 of a method's count table."""
    with open(Path(quant_dir) / TABLES[method][0], 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def common_table(methods, quant_dir, name_to_id, prefix, mirflowz_gff=None, unitas_hairpins=None, mature=None):
    unknown = [m for m in methods if m not in ADAPTERS]
    if unknown:
        raise ValueError('no benchmark adapter for method(s) {}; available: {}'.format(', '.join(unknown), ', '.join(ADAPTERS)))
    rows = []
    for method in methods:
        rows += ADAPTERS[method](quant_dir, name_to_id, prefix=prefix, mirflowz_gff=mirflowz_gff, unitas_hairpins=unitas_hairpins, mature=mature)
    return [row for row in rows if row[4]]
