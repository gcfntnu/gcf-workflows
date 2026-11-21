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
- biocore_htseq (the genome counter of an earlier in-house pipeline, `htseq_counts.tsv`, smallrna/rules/quant/biocore_htseq.smk)
  reports miRNAs by miRBase name, isomiRs summed; a read counts in full for every name it overlaps.
- biocore_miraligner (`mirna.tsv`, smallrna/rules/quant/biocore_miraligner.smk) reports miRNAs by miRBase name, its substitution-free
  isomiRs summed; a sequence annotated to several mature miRNAs is divided between them.

Adding a method means adding one function to ADAPTERS that reads its output directory, and its count table to TABLES:
before scoring, benchmark.py checks that table's samples and the miRBase release in the method's reference.tsv (#407).

isomiR table (#454): for the methods that report read sequences, one row per method, sample, reported miRNA and read
sequence (ISOMIR_COLUMNS), counted by the method's own rule, so that a miRNA's rows add up to its count in the common
table (check_isomirs). `label` is the method's own call, canonical or isomiR (empty for UNITAS, which has none):

- UNITAS: isomir_counts.tsv (`sequence|name`, unitas_isomirtable.py), the name without the modification UNITAS appends
  to a read with a substitution (`miR-548x-3p(1T->A)`), as mir_counts.tsv counts it; names placed as in the common table;
- miRge3: mapped.csv, as rules/quant/scripts/mirge3_counts.py counts miR.Counts.csv: names folded by the project
  library's merge file, a merge group's reads on its first member, and only miRNAs with at least MIRGE3_MIN_EXACT exact
  reads in the project; a read under `exact miRNA` is canonical;
- mirflowz: the `name|5' shift|3' shift|CIGAR|MD|sequence` rows are isomiRs; a plain row counts the reads that match
  mirflowz's annotated mature exactly, whose sequence the table does not give (empty);
- nf-core/smrnaseq: joined_samples_mirtop.tsv, the first row of each UID (rules/quant/scripts/smrnaseq_counts.py);
  `Variant` NA is canonical;
- biocore_miraligner: each sample's isomir.tsv (`<miRNA>_<sequence>_<mism;add;t5;t3>`, a sequence's reads already
  divided by its ambiguity), the rows without a substitution, as its mature table counts them; canonical when it has no
  addition and no 5' or 3' trimming.

TEsmall, tiny-count and biocore_htseq report features only.
"""
import collections
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


def biocore_htseq(quant_dir, name_to_id, **_):
    """The biocore_htseq route's joined mature counts (miRNA, samples), by miRBase name."""
    return named_counts('biocore_htseq', Path(quant_dir) / 'biocore_htseq' / 'htseq_counts.tsv', name_to_id)


def biocore_miraligner(quant_dir, name_to_id, **_):
    """The biocore_miraligner route's joined mature counts (miRNA, samples), by miRBase name."""
    return named_counts('biocore_miraligner', Path(quant_dir) / 'biocore_miraligner' / 'mirna.tsv', name_to_id)


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
            'biocore_htseq': biocore_htseq, 'biocore_miraligner': biocore_miraligner}

# Each method's count table under the project's quant directory: (path, delimiter, columns before the samples).
TABLES = {'unitas': ('unitas/mir_counts.tsv', '\t', 1), 'mirge3': ('mirge3/miR.Counts.csv', ',', 1),
          'mirflowz': ('mirflowz/all_mirna_counts.tab', '\t', 1), 'tesmall': ('tesmall/tesmall_counts.tsv', '\t', 2),
          'tinyrna': ('tinyrna/tinyrna_counts.tsv', '\t', 3), 'smrnaseq': ('smrnaseq/mirna.tsv', '\t', 1),
          'biocore_htseq': ('biocore_htseq/htseq_counts.tsv', '\t', 1), 'biocore_miraligner': ('biocore_miraligner/mirna.tsv', '\t', 1)}
# Methods whose reference.tsv is written by quant_reference after their count table, with that table's SHA-256.
TABLE_BOUND_RECORDS = ('mirflowz', 'tesmall', 'tinyrna', 'smrnaseq', 'biocore_htseq', 'biocore_miraligner')
# Methods renamed since (#434). A project prepared earlier names them in its benchmark/dataset.yaml; evaluating it stops
# with the new names rather than reading the outputs under the old ones.
RENAMED = {'htseq': 'biocore_htseq', 'miraligner': 'biocore_miraligner'}


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


# isomiR table (#454)
ISOMIR_COLUMNS = ('method', 'sample', 'mature_name', 'mature_id', 'sequence', 'label', 'count')
CANONICAL, ISOMIR = 'canonical', 'isomiR'
UNITAS_MODIFICATION = re.compile(r'(.*)\(.*\)\s*$')
MIRGE3_DATABASES = {'mirbase': 'miRBase', 'mirgenedb': 'MirGeneDB'}  # the merge files' names (rules/quant/mirge3.smk)
MIRGE3_MIN_EXACT = 2  # MIN_PROJECT_EXACT of rules/quant/scripts/mirge3_counts.py
MIRGE3_ANNOTATION = {'Sequence', 'annotFlag', 'exact miRNA', 'hairpin miRNA', 'mature tRNA', 'primary tRNA', 'snoRNA', 'rRNA',
                     'ncrna others', 'mRNA', 'isomiR miRNA', 'spike-in'}  # mapped.csv's columns other than the samples
MIRTOP_ANNOTATION = ('UID', 'Read', 'miRNA', 'Variant', 'iso_5p', 'iso_3p', 'iso_add3p', 'iso_snp', 'iso_5p_nt', 'iso_3p_nt',
                     'iso_add3p_nt', 'iso_snp_nt')


def unitas_isomirs(quant_dir, name_to_id, prefix, unitas_hairpins=None, mature=None, **_):
    rows, placed = [], {}
    for feature, counts in read_matrix(Path(quant_dir) / 'unitas' / 'isomir_counts.tsv', '\t').items():
        sequence, name = feature.split('|', 1)
        modified = UNITAS_MODIFICATION.match(name)
        name = (modified.group(1) if modified else name).strip()
        if name not in placed:
            placed[name] = unitas_place(name, name_to_id, prefix, unitas_hairpins or {}, mature or {})[0]
        rows += [('unitas', s, name, placed[name], sequence, '', c) for s, c in counts.items()]
    return rows


def mirge3_isomirs(quant_dir, name_to_id, **_):
    base = Path(quant_dir) / 'mirge3'
    database = method_record('mirge3', quant_dir, ('mirna_database',))['mirna_database']
    merge_files = sorted(base.glob('lib/*/annotation.Libs/*_merges_{}.csv'.format(MIRGE3_DATABASES.get(database.lower(), database))))
    if len(merge_files) != 1:
        raise ValueError('{}: expected one merge file of the {} library under lib/*/annotation.Libs, found {}'.format(
            base, database, len(merge_files)))
    groups = {}
    with open(merge_files[0]) as fh:
        for line in fh:
            fields = line.strip().split(',')
            groups.update((member, fields[0]) for member in fields[1:])
    rows, exact = [], collections.Counter()
    with open(base / 'mapped.csv', newline='') as fh:
        reader = csv.DictReader(fh)
        samples = [c for c in reader.fieldnames or [] if c not in MIRGE3_ANNOTATION]
        for row in reader:
            exact_name, isomir_name = row['exact miRNA'], row['isomiR miRNA']
            if exact_name and isomir_name:
                raise ValueError('{}: read {} has both an exact ({}) and an isomiR ({}) miRNA'.format(
                    base / 'mapped.csv', row['Sequence'], exact_name, isomir_name))
            if not (exact_name or isomir_name):
                continue
            name = groups.get(exact_name or isomir_name, exact_name or isomir_name)
            for sample in samples:
                count = float(row[sample])
                if exact_name:
                    exact[name] += count
                if count:
                    rows.append(('mirge3', sample, name, name_to_id.get(name.split('/')[0], ''), row['Sequence'],
                                 CANONICAL if exact_name else ISOMIR, count))
    return [row for row in rows if exact[row[2]] >= MIRGE3_MIN_EXACT]


def mirflowz_isomirs(quant_dir, name_to_id, mirflowz_gff=None, **_):
    locus_names = mirflowz_names(mirflowz_gff) if mirflowz_gff else {}
    rows = []
    for feature, counts in read_matrix(Path(quant_dir) / 'mirflowz' / 'all_mirna_counts.tab', '\t').items():
        fields = feature.split('|')
        name = fields[0]
        sequence, label = ('', CANONICAL) if len(fields) == 1 else (fields[-1], ISOMIR)
        accession = locus_names.get(name) or name_to_id.get(name, '')
        rows += [('mirflowz', s, name, accession, sequence, label, c) for s, c in counts.items()]
    return rows


def smrnaseq_isomirs(quant_dir, name_to_id, **_):
    path = Path(quant_dir) / 'smrnaseq' / 'joined_samples_mirtop.tsv'
    rows, seen = [], set()
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, [])
        if tuple(header[:len(MIRTOP_ANNOTATION)]) != MIRTOP_ANNOTATION or len(header) == len(MIRTOP_ANNOTATION):
            raise ValueError('{}: expected the columns {} and the samples, got {}'.format(path, ', '.join(MIRTOP_ANNOTATION), header))
        samples = header[len(MIRTOP_ANNOTATION):]
        for i, row in enumerate(reader, 2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError('{}:{}: {} fields, expected {}'.format(path, i, len(row), len(header)))
            if row[0] in seen:  # one collapsed read annotated twice: the route counts its first row only
                continue
            seen.add(row[0])
            label = CANONICAL if row[3] == 'NA' else ISOMIR
            rows += [('smrnaseq', s, row[2], name_to_id.get(row[2], ''), row[1], label, float(c))
                     for s, c in zip(samples, row[len(MIRTOP_ANNOTATION):])]
    return rows


def biocore_miraligner_isomirs(quant_dir, name_to_id, **_):
    rows = []
    for sample in table_samples('biocore_miraligner', quant_dir):
        path = Path(quant_dir) / 'biocore_miraligner' / sample / 'isomir.tsv'
        totals = collections.defaultdict(float)
        with open(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            if header != ['isomiR', 'count']:
                raise ValueError('{}: expected the columns isomiR, count, got {}'.format(path, header))
            for i, line in enumerate(fh, 2):
                fields = line.rstrip('\n').split('\t')
                parts = fields[0].rsplit('_', 2)
                variant = parts[-1].split(';')
                if len(fields) != 2 or len(parts) != 3 or len(variant) != 4:
                    raise ValueError('{}:{}: expected `<miRNA>_<sequence>_<mism;add;t5;t3>` and a count, got {!r}'.format(path, i, line[:100]))
                if variant[0] != '0':  # a substitution: the route leaves it out of its mature table
                    continue
                label = CANONICAL if variant[1:] == ['0', '0', '0'] else ISOMIR
                totals[(parts[0], parts[1], label)] += float(fields[1])
        rows += [('biocore_miraligner', sample, name, name_to_id.get(name, ''), sequence, label, count)
                 for (name, sequence, label), count in totals.items()]
    return rows


ISOMIR_ADAPTERS = {'unitas': unitas_isomirs, 'mirge3': mirge3_isomirs, 'mirflowz': mirflowz_isomirs, 'smrnaseq': smrnaseq_isomirs,
                   'biocore_miraligner': biocore_miraligner_isomirs}


def isomir_table(methods, quant_dir, name_to_id, prefix, mirflowz_gff=None, unitas_hairpins=None, mature=None):
    """(rows of ISOMIR_COLUMNS, methods without isomiR output): the rows of each method in ISOMIR_ADAPTERS, summed per
    method, sample, name, accession, sequence and label, without zero counts."""
    totals = collections.defaultdict(float)
    for method in methods:
        if method in ISOMIR_ADAPTERS:
            for *key, count in ISOMIR_ADAPTERS[method](quant_dir, name_to_id, prefix=prefix, mirflowz_gff=mirflowz_gff,
                                                        unitas_hairpins=unitas_hairpins, mature=mature):
                totals[tuple(key)] += count
    return [key + (count,) for key, count in totals.items() if count], [m for m in methods if m not in ISOMIR_ADAPTERS]


def check_isomirs(isomirs, common, methods):
    """ValueError unless, for each of `methods`, the isomiR rows of each sample and reported miRNA add up to its count
    in the common table (relative tolerance 1e-6): the isomiR table must be the method's own counting."""
    sums, counts = collections.defaultdict(float), collections.defaultdict(float)
    for method, sample, name, _, _, _, count in isomirs:
        sums[(method, sample, name)] += count
    for method, sample, name, _, count in common:
        if method in methods:
            counts[(method, sample, name)] += float(count)
    wrong = sorted(key for key in set(sums) | set(counts) if not math.isclose(sums[key], counts[key], rel_tol=1e-6, abs_tol=1e-6))
    if wrong:
        method, sample, name = wrong[0]
        raise ValueError('{} isomiR rows do not add up to the common table ({} of {} method, sample and miRNA), e.g. {} {} {}: {:g} in '
                         'the isomiR rows, {:g} in the common table'.format(', '.join(sorted({w[0] for w in wrong})), len(wrong),
                                                                           len(set(sums) | set(counts)), method, sample, name,
                                                                           sums[wrong[0]], counts[wrong[0]]))
