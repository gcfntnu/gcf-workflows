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
  the same mature miRNA (`MIMAT0000416_1`); the loci are summed per accession.

Adding a method means adding one function to ADAPTERS that reads its output directory.
"""
import csv
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


ADAPTERS = {'unitas': unitas, 'mirge3': mirge3, 'mirflowz': mirflowz, 'tesmall': tesmall, 'tinyrna': tinyrna}


def common_table(methods, quant_dir, name_to_id, prefix, mirflowz_gff=None, unitas_hairpins=None, mature=None):
    unknown = [m for m in methods if m not in ADAPTERS]
    if unknown:
        raise ValueError('no benchmark adapter for method(s) {}; available: {}'.format(', '.join(unknown), ', '.join(ADAPTERS)))
    rows = []
    for method in methods:
        rows += ADAPTERS[method](quant_dir, name_to_id, prefix=prefix, mirflowz_gff=mirflowz_gff, unitas_hairpins=unitas_hairpins, mature=mature)
    return [row for row in rows if row[4]]
