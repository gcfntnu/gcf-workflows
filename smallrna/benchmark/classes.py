"""RNA classes and read lengths across routes (#504).

Each route's own class assignment, read from the files it writes per read sequence, mapped to one common set of
classes, with read lengths, so routes can be compared beyond miRNA. Nothing is quantified again.

Routes and their sources (`quant/<route>/`):
- UNITAS: `unitas/<sample>/unitas.full_annotation_matrix.txt` (sequence, reads, annotations). A sequence's reads are
  divided evenly between its annotations, which reproduces UNITAS's own `unitas.annotation_summary.txt`. An annotation
  starting with a class name (`lncRNA|…`) is that class; miRNA names are miRNA, of the organism when they carry no
  species prefix (`miR-21-5p`) and of other species when they do (`mmu-mir-…`, `osa-MIR168a-5p`), as UNITAS's
  `miRNA:<organism>` and `miRNA:other`; tRF names (`tRF-1 from …`, `5p-tR-half from …`) are tRNA, mitochondrial when
  from `MT-`; `Mt_rRNA|…` is rRNA; `piR-cluster:…` is piRNA.
- miRge3: `mapped.csv`, one class column filled per sequence (miRge3 assigns sequentially); exact and isomiR miRNA are
  miRNA. `mapped.csv` and `unmapped.csv` together hold every filtered read, so they are also the reads' lengths when
  miRge3 is evaluated (the baseline below).
- TEsmall: `<sample>/<sample>.anno.rlen.info` (its classes by read length); `<sample>.rRNA.bam`, the reads TEsmall's
  rRNA filter removed before mapping (rRNA, by length); `count_summary.txt` for its structural RNA by type (`tRNA:`,
  `snoRNA`, …), totals only: TEsmall's lengths have structural RNA as one class. Its 3' tRF pass is left out.
- tiny-count: `<sample>/assigned_nt_len_dist.csv`, its assigned reads (miRNA features only) by length.
- nf-core/smrnaseq: `joined_samples_mirtop.tsv`, or each sample's `mirtop/<sample>/counts/mirtop.tsv` where a
  project's cleanup removed the joined table; one row per collapsed read, its first, as the route counts (adapters.py).
- biocore_miraligner: its isomiR tables (adapters.py), a read divided between its miRNAs; miRNA only.
- biocore_htseq: `<sample>/isomir.tsv` (miRBase miRNAs) and, where the route writes them (#503), `<sample>/ncrna.tsv`,
  each sequence's reads divided between the distinct types of its features, as the earlier pipeline's class totals.
  A read can be in both (the pipeline's split per alignment record), so its unassigned reads are those of sequences
  in neither table.

Common classes (CLASSES) and states per route (STATES): `reported`; `lumped` into a coarser class the route reports
instead (miRge3's `ncrna others`, TEsmall's structural RNA and gene exons); `totals only` (TEsmall's structural RNA by
type); `not reported`; `not available` (a route's file missing, e.g. TEsmall's rRNA BAM). A native class the mapping does not know becomes `unmapped`, with its name kept.

Baseline: each sample's filtered reads by length, from miRge3's tables or else the filtered FASTQ; `unassigned` is the
baseline minus the reads a route assigned. Shares are of the baseline's reads.
"""
import collections
import csv
import gzip
import re
import struct
from pathlib import Path

import adapters

CLASSES = ('miRNA', 'miRNA (other species)', 'hairpin', 'tRNA', 'rRNA', 'snoRNA', 'snRNA', 'misc ncRNA', 'piRNA',
           'lncRNA', 'mRNA', 'intron', 'pseudogene', 'TE', 'low complexity', 'ncRNA (miRge3, lumped)',
           'structural RNA (TEsmall, lumped)', 'gene exons (TEsmall, lumped)', 'other', 'unmapped', 'unassigned')
LUMP_MIRGE3, LUMP_TESMALL, LUMP_EXONS = 'ncRNA (miRge3, lumped)', 'structural RNA (TEsmall, lumped)', 'gene exons (TEsmall, lumped)'
ROUTES = ('unitas', 'mirge3', 'tesmall', 'tinyrna', 'smrnaseq', 'biocore_htseq', 'biocore_miraligner')

NATIVE = {
    'miRNA': 'miRNA', 'exact miRNA': 'miRNA', 'isomiR miRNA': 'miRNA', 'miRBase miRNA': 'miRNA',
    'miRNA (other species)': 'miRNA (other species)',
    'hairpin': 'hairpin', 'hairpin miRNA': 'hairpin', 'pre_miRNA': 'hairpin', 'precursor_RNA': 'hairpin',
    'miRNA (RNAcentral/Ensembl)': 'hairpin',
    'tRNA': 'tRNA', 'Mt_tRNA': 'tRNA', 'mature tRNA': 'tRNA', 'primary tRNA': 'tRNA',
    'rRNA': 'rRNA', 'Mt_rRNA': 'rRNA', 'rRNA_pseudogene': 'rRNA', 'rRNA (pre-filter)': 'rRNA',
    'snoRNA': 'snoRNA', 'scaRNA': 'snoRNA', 'snRNA': 'snRNA',
    'misc_RNA': 'misc ncRNA', 'Y_RNA': 'misc ncRNA', 'vault_RNA': 'misc ncRNA', 'scRNA': 'misc ncRNA',
    'SRP_RNA': 'misc ncRNA', 'RNase_P_RNA': 'misc ncRNA', 'RNase_MRP_RNA': 'misc ncRNA', 'ribozyme': 'misc ncRNA',
    'hammerhead_ribozyme': 'misc ncRNA', 'telomerase_RNA': 'misc ncRNA', 'sRNA': 'misc ncRNA', 'guide_RNA': 'misc ncRNA',
    'autocatalytically_spliced_intron': 'misc ncRNA', '7SK': 'misc ncRNA', '7SK_RNA': 'misc ncRNA',
    'piRNA': 'piRNA', 'piRNA_cluster': 'piRNA', 'piRNA cluster': 'piRNA',
    'lncRNA': 'lncRNA', 'lincRNA': 'lncRNA', 'antisense': 'lncRNA', 'antisense_RNA': 'lncRNA',
    'processed_transcript': 'lncRNA', 'sense_intronic': 'lncRNA', 'sense_overlapping': 'lncRNA',
    '3prime_overlapping_ncRNA': 'lncRNA', 'bidirectional_promoter_lncRNA': 'lncRNA', 'macro_lncRNA': 'lncRNA',
    'non_coding': 'lncRNA', 'TEC': 'lncRNA',
    'protein_coding': 'mRNA', 'mRNA': 'mRNA', 'exon': LUMP_EXONS, 'nonsense_mediated_decay': 'mRNA',
    'non_stop_decay': 'mRNA', 'retained_intron': 'mRNA', 'protein_coding_CDS_not_defined': 'mRNA',
    'protein_coding_LoF': 'mRNA',
    'intron': 'intron', 'sense_TE': 'TE', 'anti_TE': 'TE', 'low_complexity': 'low complexity',
    'ncrna others': LUMP_MIRGE3, 'structural_RNA': LUMP_TESMALL,
    'ncRNA': 'other', 'other': 'other', 'spike-in': 'other', 'artifact': 'other', 'misc_rna': 'other',
}
MIRNA_WORD = r'(miR|mir|MIR|let|lin|lsy|bantam|iab)(-|(?=[0-9]))'  # plant names: `MIR156a`
SPECIES_PREFIX = re.compile(r'^[a-z][a-z0-9]{1,5}-' + MIRNA_WORD)  # `mmu-mir-…`, `hsv1-mir-…`, `osa-MIR…`; not `let-7a-5p`
MIRNA_NAME = re.compile(r'^([a-z][a-z0-9]{1,5}-)?' + MIRNA_WORD)
TRF_NAME = re.compile(r"^(\d+p-tR-half|tRF-1|misc-tRF|tRNA-leader|\d'tRF|\d'CCA-tRF|\d'tR-half|\d+p-tRF|3p-CCA-tRF)")

ALL = set(CLASSES)
STATES = {
    'unitas': {c: 'reported' for c in ('miRNA', 'miRNA (other species)', 'tRNA', 'rRNA', 'snoRNA', 'snRNA', 'misc ncRNA',
                                        'piRNA', 'lncRNA', 'mRNA', 'pseudogene', 'low complexity', 'other')},
    'mirge3': {**{c: 'reported' for c in ('miRNA', 'hairpin', 'tRNA', 'rRNA', 'snoRNA', 'mRNA', LUMP_MIRGE3, 'other')},
               **{c: 'lumped in ' + LUMP_MIRGE3 for c in ('snRNA', 'misc ncRNA', 'lncRNA', 'pseudogene')}},
    'tesmall': {**{c: 'reported' for c in ('miRNA', 'hairpin', 'rRNA', 'piRNA', 'intron', 'TE', LUMP_TESMALL, LUMP_EXONS)},
                **{c: 'totals only (in {})'.format(LUMP_TESMALL) for c in ('tRNA', 'snoRNA', 'snRNA', 'misc ncRNA')},
                **{c: 'lumped in ' + LUMP_EXONS for c in ('mRNA', 'lncRNA', 'pseudogene')}},
    'tinyrna': {'miRNA': 'reported'}, 'smrnaseq': {'miRNA': 'reported'}, 'biocore_miraligner': {'miRNA': 'reported'},
    'biocore_htseq': {c: 'reported' for c in ('miRNA', 'hairpin', 'tRNA', 'rRNA', 'snoRNA', 'snRNA', 'misc ncRNA', 'lncRNA',
                                               'mRNA', 'pseudogene', 'other')},
}


def common(native):
    """The common class of a native class name; `unmapped` for one the mapping does not know."""
    if native in NATIVE:
        return NATIVE[native]
    if native.endswith('pseudogene'):
        return 'pseudogene'
    if re.match(r'^(IG|TR)_[A-Z]+_gene$', native):
        return 'mRNA'
    if native.startswith('structural_RNA:'):  # TEsmall's structural RNA by type (totals only)
        return common(native.split(':', 1)[1])
    return 'unmapped'


def _float(text):
    return float(text) if text not in ('', 'NA') else 0.0


class Counts:
    """A route's reads per sample: by (native class, length) (length None: a total without lengths) and the reads it
    assigned at all, by length (for `unassigned`)."""

    def __init__(self):
        self.classes = collections.defaultdict(collections.Counter)
        self.assigned = collections.defaultdict(collections.Counter)
        self.unavailable = {}  # {class: reason}: a file the route writes is missing, so its reads would look unassigned

    def add(self, sample, native, length, reads):
        if reads:
            self.classes[sample][(native, length)] += reads


# --- the routes

def unitas_native(head, vocabulary):
    if head in vocabulary and head != 'miRNA':
        return head
    if TRF_NAME.match(head) or ' from tRNA-' in head or ' from MT-' in head or ' from mt-' in head:
        return 'Mt_tRNA' if (' from MT-' in head or ' from mt-' in head) else 'tRNA'
    if head == 'Mt_rRNA':
        return 'Mt_rRNA'
    if head.startswith('piR-cluster') or head.startswith('piR-'):
        return 'piRNA cluster'
    if MIRNA_NAME.match(head):
        return 'miRNA (other species)' if SPECIES_PREFIX.match(head) else 'miRNA'
    return head


def unitas(quant, samples, skip=frozenset()):
    out = Counts()
    for sample in samples:
        folder = Path(quant) / 'unitas' / 'unitas' / sample
        vocabulary = set()
        with open(folder / 'unitas.annotation_summary.txt') as fh:
            for line in fh:
                if line.strip() and not line[0].isspace():
                    vocabulary.add(line.split('\t')[0])
        with open(folder / 'unitas.full_annotation_matrix.txt') as fh:
            next(fh, None)
            for line in fh:
                fields = line.rstrip('\n').split('\t')
                annotations = [a for a in fields[2:] if a]
                if len(fields) < 3 or not annotations or fields[0] in skip:
                    continue
                reads, length = _float(fields[1]), len(fields[0])
                out.assigned[sample][length] += reads
                for a in annotations:
                    out.add(sample, unitas_native(a.split('|')[0].strip(), vocabulary), length, reads / len(annotations))
    return out


COORDINATES = re.compile(r'\s*\(\d+-\d+\)$')       # ` (517-535)`: where the read lies in the record
ENSEMBL_GENE = re.compile(r'ENS[A-Z]*G\d+')
SILVA_SUBUNIT = re.compile(r'^(SSU|LSU):')          # SILVA records: one per genomic copy, counted by subunit
PIRNA_CLUSTER = re.compile(r'^piR-cluster:([^:]+)')  # `piR-cluster:Hsap_0125_2:24815296..24831604(17A->C)`
MISMATCH = re.compile(r'\((\d+[ACGTN]->[ACGTN],?)+\)$')
TRANSCRIPT = re.compile(r'[-_]\(\d+\)$')  # a gene's further transcripts: '<symbol>-<gene>_(n)', in dumps before 20-11-2026 '<symbol>-(n)'  # `miR-21-5p(12C->A)`: a hit with a mismatch, the same feature
UNASSIGNED = ('unassigned', 'unassigned')


def unitas_feature(annotation, vocabulary):
    """(common class, feature) of one UNITAS annotation of a read sequence (#409): the feature a class's detection
    counts. Ensembl records (`class|NAME-ENSG… (a-b)`) and tRNA fragments (`5p-tRF from MT-TS1-ENSG… (a-b)`) count
    by Ensembl gene where the record names one, else by record name without coordinates; SILVA rRNA copies by
    subunit (SSU, LSU); piRNA clusters by cluster; miRNAs and the rest by name. A mismatch UNITAS notes after a name
    (`miR-21-5p(12C->A)`) is the same feature."""
    annotation = annotation.strip()
    head, bar, record = annotation.partition('|')
    head = head.strip()
    cls = common(unitas_native(head, vocabulary))
    if bar:
        name = COORDINATES.sub('', record.strip())
        subunit = SILVA_SUBUNIT.match(name)
        if subunit:
            return cls, subunit.group(1)
    elif ' from ' in head:
        name = MISMATCH.sub('', COORDINATES.sub('', head.split(' from ', 1)[1]))
    else:
        cluster = PIRNA_CLUSTER.match(head)
        return cls, cluster.group(1) if cluster else MISMATCH.sub('', COORDINATES.sub('', head))
    gene = ENSEMBL_GENE.search(name)
    return cls, gene.group(0) if gene else TRANSCRIPT.sub('', name)


def unitas_sequences(quant, sample):
    """Iterator of (reads, {(class, feature): share}) over a sample's UNITAS read sequences (#409), from its full
    annotation matrix: each annotated sequence's reads with the share each feature gets under UNITAS's equal split over
    its annotations (one feature's records together). A sequence without annotation is ('unassigned', 'unassigned'),
    so the reads of all rows are UNITAS's input reads."""
    folder = Path(quant) / 'unitas' / 'unitas' / sample
    vocabulary = set()
    with open(folder / 'unitas.annotation_summary.txt') as fh:
        for line in fh:
            if line.strip() and not line[0].isspace():
                vocabulary.add(line.split('\t')[0])
    cache = {}  # by the annotation without its coordinates, so one entry per record rather than per read position
    with open(folder / 'unitas.full_annotation_matrix.txt') as fh:
        if 'annotation(s)' not in next(fh, ''):  # as the workflow's parser (unitas_feature_counts.py): empty or not a matrix
            raise ValueError('{}: no UNITAS full annotation matrix header'.format(folder / 'unitas.full_annotation_matrix.txt'))
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                continue
            annotations = [a for a in fields[2:] if a.strip()]
            hits = collections.Counter()
            for a in annotations:
                record = a.rsplit(' (', 1)[0] if a.endswith(')') and ' (' in a else a
                key = cache.get(record)
                if key is None:
                    key = cache[record] = unitas_feature(record, vocabulary)
                hits[key] += 1
            # one division per feature: six annotations of one feature give it 1.0, not six sixths summed (0.999...)
            yield int(round(_float(fields[1]))), {k: n / len(annotations) for k, n in hits.items()} or {UNASSIGNED: 1.0}


def mirge3_tables(quant):
    """(mapped rows, unmapped rows, class columns, sample columns) of miRge3's mapped.csv and unmapped.csv."""
    tables = []
    for name in ('mapped.csv', 'unmapped.csv'):
        path = Path(quant) / 'mirge3' / name
        try:
            with open(path, newline='') as fh:
                reader = csv.reader(fh)
                header = next(reader)
                tables.append((header, list(reader)))
        except (OSError, csv.Error, StopIteration):  # missing, empty or unreadable: each table on its own
            tables.append(None)
    return tables


def mirge3(quant, samples, tables=None, skip=frozenset()):
    out = Counts()
    mapped, _ = tables or mirge3_tables(quant)
    if mapped is None:
        raise FileNotFoundError(str(Path(quant) / 'mirge3' / 'mapped.csv'))
    header, rows = mapped
    classes = [i for i, c in enumerate(header) if c in adapters.MIRGE3_ANNOTATION and c not in ('Sequence', 'annotFlag')]
    columns = {s: header.index(s) for s in samples if s in header}
    for row in rows:
        natives = [header[i] for i in classes if row[i]]
        if not natives or row[0] in skip:
            continue
        length = len(row[0])
        for sample, i in columns.items():
            reads = _float(row[i])
            if reads:
                out.assigned[sample][length] += reads
                for native in natives:
                    out.add(sample, native, length, reads / len(natives))
    return out


def bam_read_lengths(path):
    """{length: reads} of a BAM's mapped primary records, each read name once (a BGZF file is gzip)."""
    lengths, seen = collections.Counter(), set()
    with gzip.open(path, 'rb') as fh:
        if fh.read(4) != b'BAM\x01':
            raise ValueError('{}: not a BAM file'.format(path))
        fh.read(struct.unpack('<i', fh.read(4))[0])
        for _ in range(struct.unpack('<i', fh.read(4))[0]):
            fh.read(struct.unpack('<i', fh.read(4))[0] + 4)
        while True:
            head = fh.read(4)
            if len(head) < 4:
                break
            data = fh.read(struct.unpack('<i', head)[0])
            _, _, l_name, _, _, _, flag, l_seq = struct.unpack('<iiBBHHHi', data[:20])
            if flag & 0x904:  # unmapped, secondary, supplementary
                continue
            name = data[32:32 + l_name - 1]
            if name not in seen:
                seen.add(name)
                lengths[l_seq] += 1
    return lengths


def tesmall(quant, samples):
    out = Counts()
    for sample in samples:
        folder = Path(quant) / 'tesmall' / sample
        with open(folder / '{}.anno.rlen.info'.format(sample), newline='') as fh:
            rows = list(csv.reader(fh, delimiter='\t'))
        header = rows[0]
        for row in rows[1:]:
            length = int(row[0])
            for native, value in zip(header[1:], row[1:]):
                reads = _float(value)
                out.add(sample, native, length, reads)
                out.assigned[sample][length] += reads
        rrna = folder / '{}.rRNA.bam'.format(sample)
        if rrna.is_file():
            for length, reads in bam_read_lengths(rrna).items():
                out.add(sample, 'rRNA (pre-filter)', length, reads)
                out.assigned[sample][length] += reads
        else:  # the reads its rRNA filter removed are unknown: neither rRNA nor unassigned can be given
            out.unavailable['rRNA'] = out.unavailable['unassigned'] = 'not available ({} missing)'.format(rrna.name)
        summary = folder / 'count_summary.txt'
        if summary.is_file():
            with open(summary, newline='') as fh:
                reader = csv.reader(fh, delimiter='\t')
                next(reader, None)
                for row in reader:
                    if len(row) >= 3 and row[1] == 'structural_RNA':
                        out.add(sample, 'structural_RNA:' + re.split(r'[:]', row[0])[0], None, _float(row[2]))
    return out


def tinyrna(quant, samples, skip=frozenset()):  # no sequences: `skip` cannot apply
    out = Counts()
    for sample in samples:
        with open(Path(quant) / 'tinyrna' / sample / 'assigned_nt_len_dist.csv', newline='') as fh:
            reader = csv.reader(fh)
            next(reader, None)
            for row in reader:
                length, reads = int(row[0]), sum(_float(v) for v in row[1:])
                out.add(sample, 'miRNA', length, reads)
                out.assigned[sample][length] += reads
    return out


def isomir_route(method):
    def parse(quant, samples, skip=frozenset()):
        out = Counts()
        per_sequence = collections.defaultdict(float)
        for _, sample, _, _, sequence, _, count in adapters.ISOMIR_ADAPTERS[method](quant, {}):
            per_sequence[(sample, sequence)] += count  # smrnaseq: one row per read; miraligner: a read divided between names
        for (sample, sequence), reads in per_sequence.items():
            if sample in samples and sequence not in skip:
                out.add(sample, 'miRNA', len(sequence), reads)
                out.assigned[sample][len(sequence)] += reads
        return out
    return parse


def keyed_table(path):
    """{sequence: (names, count)} of a `<name>_<sequence>`, count table; a sequence's count is the same for each of its
    names (each (name, read) pair counts once), so the largest is kept."""
    out = {}
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            key, count = line.rstrip('\n').rsplit('\t', 1)
            name, sequence = key.rsplit('_', 1)
            names, value = out.get(sequence, (set(), 0.0))
            names.add(name)
            out[sequence] = (names, max(value, float(count)))
    return out


def biocore_htseq(quant, samples, skip=frozenset()):
    out = Counts()
    for sample in samples:
        folder = Path(quant) / 'biocore_htseq' / sample
        assigned = {}
        for sequence, (_, reads) in keyed_table(folder / 'isomir.tsv').items():
            if sequence in skip:
                continue
            out.add(sample, 'miRBase miRNA', len(sequence), reads)
            assigned[sequence] = reads
        if (folder / 'ncrna.tsv').is_file():
            for sequence, (names, reads) in keyed_table(folder / 'ncrna.tsv').items():
                if sequence in skip:
                    continue
                types = sorted({name.split('-', 1)[1] if '-' in name else name for name in names})
                for kind in types:
                    out.add(sample, 'miRNA (RNAcentral/Ensembl)' if kind == 'miRNA' else kind, len(sequence), reads / len(types))
                assigned[sequence] = max(reads, assigned.get(sequence, 0.0))
        for sequence, reads in assigned.items():  # the union: a read in both paths is assigned once
            out.assigned[sample][len(sequence)] += reads
    return out


def smrnaseq(quant, samples, skip=frozenset()):
    out = Counts()
    joined = Path(quant) / 'smrnaseq' / 'joined_samples_mirtop.tsv'
    paths = [joined] if joined.is_file() else [Path(quant) / 'smrnaseq' / 'mirtop' / s / 'counts' / 'mirtop.tsv' for s in samples]
    for path in paths:
        seen = set()
        with open(path, newline='') as fh:
            reader = csv.reader(fh, delimiter='\t')
            header = next(reader, [])
            if tuple(header[:len(adapters.MIRTOP_ANNOTATION)]) != adapters.MIRTOP_ANNOTATION:
                raise ValueError('{}: expected the mirtop columns, got {}'.format(path, header[:4]))
            columns = {s: i for i, s in enumerate(header) if i >= len(adapters.MIRTOP_ANNOTATION) and s in samples}
            for row in reader:
                if not row or row[0] in seen:  # one collapsed read annotated twice: the route counts its first row only
                    continue
                seen.add(row[0])
                if row[1] in skip:
                    continue
                for sample, i in columns.items():
                    reads = _float(row[i])
                    if reads:
                        out.add(sample, 'miRNA', len(row[1]), reads)
                        out.assigned[sample][len(row[1])] += reads
    return out


PARSERS = {'unitas': unitas, 'tesmall': tesmall, 'tinyrna': tinyrna, 'smrnaseq': smrnaseq,
           'biocore_htseq': biocore_htseq, 'biocore_miraligner': isomir_route('biocore_miraligner')}


# --- the baseline and the summary

def baseline(quant, samples, tables=None, fastq=None):
    """({sample: {length: reads}} of the filtered reads, source), from miRge3's mapped.csv and unmapped.csv when both
    exist, else from the filtered FASTQ (`fastq`: {sample: path}); (None, reason) when neither is there."""
    mapped, unmapped = tables or mirge3_tables(quant)
    if mapped and unmapped:
        try:
            out = {s: collections.Counter() for s in samples}
            for header, rows in (mapped, unmapped):
                columns = {s: header.index(s) for s in samples if s in header}
                for row in rows:
                    for sample, i in columns.items():
                        reads = _float(row[i])
                        if reads:
                            out[sample][len(row[0])] += reads
            return out, "miRge3's mapped.csv and unmapped.csv"
        except (ValueError, IndexError):  # a malformed table: the filtered FASTQ instead
            pass
    if fastq and all(Path(fastq.get(s, '')).is_file() for s in samples):
        out = {}
        for sample in samples:
            lengths = collections.Counter()
            with open(fastq[sample]) as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        lengths[len(line) - 1] += 1
            out[sample] = lengths
        return out, 'the filtered FASTQ'
    return None, "neither miRge3's tables nor the filtered reads are there"


def species_mirnas_only(quant):
    """Whether UNITAS ran with -species_miR_only, as its counting settings record it (quant/unitas/settings.tsv, #417):
    then it searched no other species' miRNAs, and that class is not reported."""
    path = Path(quant) / 'unitas' / 'settings.tsv'
    try:
        return '-species_miR_only' in path.read_text()
    except OSError:
        return False


def summarize(quant, methods, samples, fastq=None):
    """(rows of classes.tsv, rows of class_states.tsv, {method: Counts}, baseline, notes). A route whose files cannot
    be read is left out with a note; the evaluation goes on."""
    tables = mirge3_tables(quant)  # the baseline too, also where miRge3 is not scored: its tables hold every filtered read
    counts, notes = {}, []
    for method in methods:
        if method not in ROUTES:
            continue
        try:
            counts[method] = mirge3(quant, samples, tables) if method == 'mirge3' else PARSERS[method](quant, samples)
        except (OSError, ValueError, IndexError, KeyError, StopIteration, struct.error) as e:
            notes.append('{}: not read ({}: {})'.format(method, type(e).__name__, e))
    try:
        base, source = baseline(quant, samples, tables, fastq)
    except (OSError, ValueError, IndexError) as e:
        base, source = None, 'not read ({}: {})'.format(type(e).__name__, e)
    rows, states = [], []
    for method, c in counts.items():
        for sample in samples:
            for (native, length), reads in sorted(c.classes[sample].items(), key=lambda kv: (kv[0][0], kv[0][1] or 0)):
                rows.append((method, sample, common(native), native, '' if length is None else length, reads))
            if base is not None and 'unassigned' not in c.unavailable:
                for length, total in sorted(base[sample].items()):
                    left = total - c.assigned[sample].get(length, 0.0)
                    if left > 1e-9:
                        rows.append((method, sample, 'unassigned', 'unassigned', length, left))
        known = STATES[method]
        if method == 'biocore_htseq' and not any(n not in ('miRBase miRNA',) for s in samples for n, _ in c.classes[s]):
            known = {'miRNA': 'reported'}  # no ncRNA tables (#503: human only)
        if method == 'unitas' and species_mirnas_only(quant):
            known = {**known, 'miRNA (other species)': 'not reported (-species_miR_only)'}
        for cls in CLASSES:
            if cls in c.unavailable:
                state = c.unavailable[cls]
            elif cls == 'unassigned':
                state = 'reported' if base is not None else 'not available (no baseline)'
            elif cls == 'unmapped':
                state = 'reported' if any(common(n) == 'unmapped' for s in samples for n, _ in c.classes[s]) else 'not reported'
            else:
                state = known.get(cls, 'not reported')
            states.append((method, cls, state))
    notes.append('baseline: ' + source)
    return rows, states, counts, base, notes


# The synthetic set's decoy probe samples (#412): each holds the carriers of P_control and the reads of one class and
# source; the common classes that are right for it (decoys.py: other_ncRNA is snRNA, 7SK, 7SL and vault RNA).
DECOY_CLASSES = {'tRNA-gtrnadb': ('tRNA',), 'tRNA-ensembl': ('tRNA',), 'rRNA-refseq': ('rRNA',), 'rRNA-ensembl': ('rRNA',),
                 'Y_RNA-ensembl': ('misc ncRNA',), 'snoRNA-ensembl': ('snoRNA',), 'other_ncRNA-ensembl': ('snRNA', 'misc ncRNA'),
                 'mRNA-ensembl': ('mRNA',), 'piRNA-pirnaclusterdb': ('piRNA',), 'precursor-mirbase': ('hairpin',)}
# Sources scored also by molecule or fragment type, as such reads in real libraries are counted: a part's pattern is
# searched in the truth's `origin` (rRNA: `<region or gene>:...`; tRNA: `<tRNA>:<fragment type>`). The decoys' mix of
# parts differs from real libraries (e.g. the 45S spacers and tRNA halves are far more common among them).
DECOY_PARTS = {'rRNA-refseq': {'18S': r'^18S:', '5.8S': r'^5\.8S:', '28S': r'^28S:', '5S': r'^5S:',
                               'spacers': r'^(5ETS|ITS1|ITS2|3ETS):'},
               'rRNA-ensembl': {'5S': r'^(RNA5S\d*|RNA5SP\d+|5S_rRNA):', '5.8S': r'^(RNA5-8SN?\d*|RNA5-8SP\d+|5_8S_rRNA):',
                                'mitochondrial': r'^MT-RNR[12]:'},
               **{source: {"5' tRF": r':5p_tRF$', "5' half": r':5p_half$', "3' tRF": r':3p_tRF$', "3' half": r':3p_half$',
                           'i-tRF': r':i_tRF$', 'modified (m1G9, m1A58)': r':(m1G9|m1A58)$'}
                  for source in ('tRNA-gtrnadb', 'tRNA-ensembl')}}


def decoy_part(source, origin):
    """The DECOY_PARTS part of a truth row's `origin`, or None."""
    for part, pattern in DECOY_PARTS.get(source, {}).items():
        if re.search(pattern, origin or ''):
            return part
    return None


DECOY_SAMPLE = re.compile(r'^P_decoy_(.+)_\d+$')
CONTROL = 'P_control'
OUTCOMES = ('right class', 'lumped', 'miRNA', 'other class', 'unassigned')


def sample_classes(rows, method, sample):
    """{class: reads} of one route and sample; TEsmall's structural RNA lump divided by its totals by type."""
    by, types = collections.Counter(), collections.Counter()
    for m, s, cls, _, length, reads in rows:
        if m == method and s == sample:
            (types if length == '' else by)[cls] += reads
    lump = by.pop(LUMP_TESMALL, 0)
    if lump and sum(types.values()):
        for cls, reads in types.items():
            by[cls] += reads * lump / sum(types.values())
    elif lump:
        by[LUMP_TESMALL] = lump
    return by


def outcome(cls, right, lumps):
    if cls in right:
        return 'right class'
    if cls in lumps:
        return 'lumped'
    if cls in ('miRNA', 'miRNA (other species)'):
        return 'miRNA'
    return 'unassigned' if cls == 'unassigned' else 'other class'


SEQUENCE_ROUTES = ('unitas', 'mirge3', 'smrnaseq', 'biocore_htseq', 'biocore_miraligner')  # outputs per read sequence


def decoy_scores(quant, rows, states, base, counts, methods, samples, carriers, parts=None):
    """([(method, source, class, reads, %)], [(method, source, outcome, %)]) of the decoy probe samples: where each route
    puts a decoy class's reads, in % of the decoy reads that passed the filter (a decoy sample's baseline minus
    P_control's). Routes with outputs per sequence are read again without the `carriers` (the truth's carrier
    sequences), since a route may place a carrier differently from sample to sample (biocore_htseq picks one of a
    read's loci); TEsmall and tiny-count, without sequences, by their classes minus P_control's. Unassigned: the decoy
    reads minus those the route assigned. `outcome`: the right class, a route's coarser class holding it (`lumped`),
    miRNA, another class, unassigned. `parts`: {source: {part: {sample: {sequence: reads}}}} of the truth (decoy_part),
    scored as `<source> <part>` for the routes with sequences, each read without the other parts' sequences and with
    the part's reads as the denominator (they all pass the filter). Empty without decoy samples or baseline."""
    decoys = collections.defaultdict(list)
    for sample in samples:
        match = DECOY_SAMPLE.match(sample)
        if match and match.group(1) in DECOY_CLASSES:
            decoys[match.group(1)].append(sample)
    probes = [CONTROL, *sum(decoys.values(), [])]
    if base is None or CONTROL not in samples or not decoys or any(s not in base for s in probes):
        return [], []
    state = {(m, c): s for m, c, s in states}
    total = {s: sum(base[s].values()) - sum(base[CONTROL].values()) for s in probes}
    detail, summary = [], []
    for method in methods:
        if method not in counts:
            continue
        if method not in SEQUENCE_ROUTES and counts[method].unavailable:  # e.g. TEsmall's rRNA BAM missing: not scored
            summary += [(method, source, name, 'n/a') for source in sorted(decoys) for name in OUTCOMES]
            continue
        if method in SEQUENCE_ROUTES:
            parsed = (mirge3 if method == 'mirge3' else PARSERS[method])(quant, probes, skip=frozenset(carriers))
            got = {s: collections.Counter() for s in probes}
            for s in probes:
                for (native, _), reads in parsed.classes[s].items():
                    got[s][common(native)] += reads
            assigned = {s: sum(parsed.assigned[s].values()) for s in probes}
        else:
            control = sample_classes(rows, method, CONTROL)
            control.pop('unassigned', None)
            got, assigned = {}, {}
            for s in probes:
                mine = sample_classes(rows, method, s)
                mine.pop('unassigned', None)
                got[s] = collections.Counter({c: mine.get(c, 0.0) - control.get(c, 0.0) for c in set(mine) | set(control)})
                assigned[s] = sum(counts[method].assigned[s].values()) - sum(counts[method].assigned[CONTROL].values())
        jobs = [(source, replicates, got, assigned, {s: total[s] for s in replicates}) for source, replicates in sorted(decoys.items())]
        if method in SEQUENCE_ROUTES:
            for source, by_part in sorted((parts or {}).items()):
                for part, reads_of in by_part.items():
                    replicates = [s for s in decoys.get(source, []) if s in reads_of]
                    others = {q for other, r in by_part.items() if other != part for s in replicates for q in r.get(s, {})}
                    parsed = (mirge3 if method == 'mirge3' else PARSERS[method])(quant, replicates, skip=frozenset(carriers) | others)
                    part_got = {s: collections.Counter() for s in replicates}
                    for s in replicates:
                        for (native, _), reads in parsed.classes[s].items():
                            part_got[s][common(native)] += reads
                    jobs.append(('{} {}'.format(source, part), replicates, part_got,
                                 {s: sum(parsed.assigned[s].values()) for s in replicates},
                                 {s: sum(reads_of[s].values()) for s in replicates}))
        for source, replicates, got, assigned, decoy_reads in jobs:
            right = DECOY_CLASSES[source.split(' ', 1)[0]]
            lumps = {state.get((method, c), '')[len('lumped in '):] for c in right
                     if state.get((method, c), '').startswith('lumped in ')}
            reads, denominator = collections.Counter(), sum(decoy_reads[s] for s in replicates)
            for s in replicates:
                reads.update(got[s])
                reads['unassigned'] += decoy_reads[s] - assigned[s]
            if denominator <= 0:
                continue
            ends = collections.Counter()
            for cls, value in sorted(reads.items()):
                if abs(value) > 1e-9:
                    detail.append((method, source, cls, round(value / len(replicates), 4), round(100 * value / denominator, 4)))
                    ends[outcome(cls, right, lumps)] += value
            reported = any(not state.get((method, c), 'not reported').startswith('not ') for c in right)
            for name in OUTCOMES:
                value = 'n/r' if name == 'right class' and not reported else round(100 * ends[name] / denominator, 2)
                summary.append((method, source, name, value))
    return detail, summary


def decoy_lines(summary, methods):
    """report.md's section on the decoy probe samples' class assignment."""
    if not summary:
        return []
    by = {(m, s, o): v for m, s, o, v in summary}
    sources = sorted({s for _, s, _ in by}, key=lambda s: (list(DECOY_CLASSES).index(s.split(' ', 1)[0]), s))
    shown = [m for m in methods if any(k[0] == m for k in by)]

    def cell(m, s, o):
        v = by.get((m, s, o))
        if v is None:
            return ''
        return v if isinstance(v, str) else '{:.1f}'.format(max(v, 0.0))
    lines = ['', '### Class assignment against the truth (#504)', '',
             "The decoy probe samples each add the reads of one class and source to P_control's carriers; a route's class "
             "reads there minus P_control's give where it puts them, in % of the decoy reads that passed the filter "
             '(mean of the replicates). Right class: one of the classes in the first column; `lumped`, the route\'s coarser '
             'class that holds it (miRge3\'s `ncrna others`, TEsmall\'s gene exons); `n/r`, the route has no such class; `n/a`, a file the route writes is missing. '
             "biocore_htseq counts a read in every class it touches, so its shares can sum to more than 100. For the routes with "
             'outputs per sequence, the rRNA decoys are also scored by molecule (RefSeq\'s 18S, 5.8S, 28S, 5S and the 45S '
             "pre-rRNA's spacers; Ensembl's 5S, 5.8S and mitochondrial rRNA) and the tRNA decoys by fragment type. The decoys' "
             'mix of these parts differs from real libraries (spacers and tRNA halves are far more common among them), and '
             'fragments beyond the annotated ends, such as the 5.8S long form\'s 5\' end, are not in the truth. Every class is '
             'in `class_truth.tsv`.', '',
             '| decoy source (right class) | ' + ' | '.join(shown) + ' |', '|' + '---|' * (len(shown) + 1)]
    for s in sources:
        cells = []
        for m in shown:
            right, lumped = cell(m, s, 'right class'), by.get((m, s, 'lumped'), 0)
            cells.append(right + (' (+{:.1f} lumped)'.format(lumped) if isinstance(lumped, float) and lumped >= 0.05 else ''))
        lines.append('| {} ({}) | {} |'.format(s, ', '.join(DECOY_CLASSES[s.split(' ', 1)[0]]), ' | '.join(cells)))
    for name, label in (('miRNA', 'counted as miRNA'), ('other class', 'put in another class'), ('unassigned', 'left unassigned')):
        lines += ['', 'Decoy reads {}, %:'.format(label), '', '| decoy source | ' + ' | '.join(shown) + ' |',
                  '|' + '---|' * (len(shown) + 1)]
        lines += ['| {} | {} |'.format(s, ' | '.join(cell(m, s, name) for m in shown)) for s in sources]
    return lines


def table_text(header, rows):
    return '\t'.join(header) + '\n' + ''.join('\t'.join(str(v) if not isinstance(v, float) else '{:.6g}'.format(v) for v in row) + '\n'
                                              for row in rows)


def shares(rows, base, methods, samples):
    """{method: {class: median over samples of the class's share of the sample's filtered reads, %}}."""
    totals = collections.defaultdict(float)
    for method, sample, cls, _, _, reads in rows:
        totals[(method, sample, cls)] += reads
    out = {}
    for method in methods:
        out[method] = {}
        for cls in CLASSES:
            values = sorted(100 * totals.get((method, s, cls), 0.0) / sum(base[s].values()) for s in samples if base and sum(base[s].values()))
            if values:
                mid = len(values) // 2
                out[method][cls] = values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2
    return out


def report_lines(rows, states, base, methods, samples, notes, plot_name=None):
    state = {(m, c): s for m, c, s in states}
    lines = ['', '## RNA classes and read lengths (#504)', '',
             "Each route's own class assignment, from the files it writes per read sequence (classes.py), mapped to common "
             'classes: the median share of each sample\'s filtered reads, %. `n/r`: the route has no such class; `lumped`: '
             'it counts these reads in a coarser class (miRge3\'s `ncrna others`, TEsmall\'s structural RNA and gene exons); `n/a`: '
             'a file the route writes is missing; '
             '`†`: TEsmall\'s structural RNA by type, totals without lengths, part of its lumped class. Classes come from each '
             "route's own references and precedence, so shares partly reflect those. Every row by sample and length is in "
             '`classes.tsv`, the states in `class_states.tsv`.', '']
    if base is None:
        return lines + ['Not computed: ' + '; '.join(notes) + '.']
    share = shares(rows, base, methods, samples)
    present = [c for c in CLASSES if any(state.get((m, c), 'not reported') != 'not reported' for m in methods)]
    lines += ['| class | ' + ' | '.join(methods) + ' |', '|' + '---|' * (len(methods) + 1)]
    for cls in present:
        cells = []
        for m in methods:
            s = state.get((m, cls), 'not reported')
            if s.startswith('not reported'):
                cells.append('n/r')
            elif s.startswith('not available'):
                cells.append('n/a')
            elif s.startswith('lumped'):
                cells.append('lumped')
            else:
                value = share.get(m, {}).get(cls)
                cells.append(('{:.2f}'.format(value) if value is not None else '0') + ('†' if s.startswith('totals only') else ''))
        lines.append('| {} | {} |'.format(cls, ' | '.join(cells)))
    lines += ['', 'Notes: ' + '; '.join(notes) + '.']
    if plot_name:
        lines += ['', 'Plots, pooled over samples as % of the filtered reads, with the classes grouped as in MultiQC\'s UNITAS section '
                  '(`other`: miRNA of other species, hairpin, piRNA, intron, pseudogene, TE, low complexity, other, unmapped): '
                  '`{}`, one bar per route by class; `{}`, reads by length per class, one line per route that reports '
                  'it.'.format(*plot_name)]
    return lines


# The plots' classes, as MultiQC's UNITAS section shows them: few groups, the rest in `other`; the detail is in the tables.
PLOT_GROUPS = (('miRNA', ('miRNA',)), ('tRNA', ('tRNA',)), ('rRNA', ('rRNA',)), ('snoRNA', ('snoRNA',)),
               ('misc ncRNA (snRNA, Y RNA, 7SK, …)', ('snRNA', 'misc ncRNA')), ('lncRNA', ('lncRNA',)), ('mRNA', ('mRNA',)),
               ('other', ('miRNA (other species)', 'hairpin', 'piRNA', 'intron', 'pseudogene', 'TE', 'low complexity', 'other',
                          'unmapped')),
               ('lumped (route-specific)', (LUMP_MIRGE3, LUMP_TESMALL, LUMP_EXONS)), ('unassigned', ('unassigned',)))
GROUP_OF = {cls: group for group, members in PLOT_GROUPS for cls in members}
SERIES = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#6250d6', '#008300', '#e34948')  # a validated categorical order
GREYS = {'other': '#7d7d78', 'lumped (route-specific)': '#b4b4ae', 'unassigned': '#dededa'}
LENGTH_PANELS = ('miRNA', 'tRNA', 'rRNA', 'snoRNA', 'misc ncRNA (snRNA, Y RNA, 7SK, …)', 'unassigned')


def plot(rows, base, methods, path, lengths_path, max_length=45):
    """Two PNGs, pooled over samples as % of the filtered reads: `path` one bar per route stacked by class group, and
    `lengths_path` reads by length per class group, one line per route (a route without the group is not drawn);
    False when matplotlib is missing."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    total = sum(sum(c.values()) for c in base.values()) or 1.0
    by_group, by_length = collections.defaultdict(float), collections.defaultdict(float)
    for method, _, cls, _, length, reads in rows:
        if length != '':  # TEsmall's totals by type are part of its lumped class
            group = GROUP_OF[common(cls) if cls not in ALL else cls]
            by_group[(method, group)] += 100 * reads / total
            by_length[(method, group, min(int(length), max_length))] += 100 * reads / total
    groups = [g for g, _ in PLOT_GROUPS if any(by_group.get((m, g)) for m in methods)]
    colour = {g: GREYS.get(g) or SERIES[i] for i, (g, _) in enumerate(PLOT_GROUPS)}
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#9a9a94',
                         'xtick.color': '#55554f', 'ytick.color': '#55554f', 'axes.labelcolor': '#33332f'})

    fig, ax = plt.subplots(figsize=(9, 0.42 * len(methods) + 1.5))
    for y, method in enumerate(methods):
        left = 0.0
        for g in groups:
            value = by_group.get((method, g), 0.0)
            if value <= 0:
                continue
            ax.barh(y, value, left=left, height=0.62, color=colour[g], edgecolor='white', linewidth=1,
                    hatch='///' if g.startswith('lumped') else None, label=g if y == 0 or g not in ax.get_legend_handles_labels()[1] else None)
            if value >= 5:
                dark = g in ('unassigned', 'lumped (route-specific)', 'tRNA', 'snoRNA', 'misc ncRNA (snRNA, Y RNA, 7SK, …)')
                ax.text(left + value / 2, y, '{:.0f}'.format(value), ha='center', va='center', fontsize=7,
                        color='#1f1f1c' if dark else 'white')
            left += value
    ax.set_yticks(range(len(methods)), methods)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel('% of filtered reads (all samples pooled)')
    ax.grid(axis='x', color='#e6e6e1', linewidth=0.6)
    ax.set_axisbelow(True)
    handles = dict(zip(*reversed(ax.get_legend_handles_labels())))
    fig.legend([handles[g] for g in groups if g in handles], [g for g in groups if g in handles], loc='upper center',
               ncol=5, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.55 / (0.42 * len(methods) + 1.5)))
    fig.savefig(path, dpi=120)
    plt.close(fig)

    panels = [g for g in LENGTH_PANELS if g in groups]
    lengths = list(range(min([16] + [k[2] for k in by_length]), max_length + 1))
    ncol = 3
    nrow = (len(panels) + ncol - 1) // ncol or 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 2.6 * nrow + 0.7), squeeze=False)
    for ax, g in zip(axes.flat, panels):
        for i, method in enumerate(methods):
            if by_group.get((method, g)):
                ax.plot(lengths, [by_length.get((method, g, n), 0.0) for n in lengths], color=SERIES[i % len(SERIES)],
                        linewidth=1.6, label=method)
        ax.set_title(g, fontsize=8.5, loc='left', color='#33332f')
        ax.grid(axis='y', color='#e6e6e1', linewidth=0.6)
        ax.set_xlim(lengths[0], max_length)
    for ax in axes.flat[len(panels):]:
        ax.axis('off')
    for ax in axes[:, 0]:
        ax.set_ylabel('% of filtered reads')
    for ax in axes[-1, :]:
        ax.set_xlabel('read length (nt; {} = {} or longer)'.format(max_length, max_length))
    fig.legend([plt.Line2D([], [], color=SERIES[i % len(SERIES)], linewidth=2) for i, _ in enumerate(methods)], methods,
               loc='upper center', ncol=len(methods), frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.45 / (2.6 * nrow + 0.7)))
    fig.savefig(lengths_path, dpi=120)
    plt.close(fig)
    return True
