#!/usr/bin/env python
"""Build the annotation and the rDNA and tDNA sequences of a TEsmall database from public, release-pinned sources.

TEsmall (`--dbfolder DIR -g NAME`) reads `DIR/genomes/NAME/`: BED6 files per class under `annotation/`, and `genome`,
`rDNA` and `tDNA` FASTA files with their `.fai` and Bowtie 1 indexes under `sequence/`. The genome FASTA and its index
are made by the gcfdb rule; this script writes the rest, on the genome's own sequence names (Ensembl), and records
where every class comes from:

- `miRNA.bed`, `hairpin.bed`: miRBase GFF features `miRNA` and `miRNA_primary_transcript`, named as miRBase does
  (`hsa-miR-21-5p`, `hsa-mir-21`), chromosomes renamed with the mirflowz chromosome map.
- `structural_RNA.bed`: Ensembl genes of the small structural biotypes (`STRUCTURAL_BIOTYPES`), named
  `<biotype>:<gene name>`, and GtRNAdb's high-confidence tRNA genes, named `tRNA:<GtRNAdb name>`.
- `exon.bed`, `intron.bed`: exons and introns (gaps between a transcript's exons) of all other Ensembl genes except
  miRNA genes, named by gene name (or ID).
- `TE.bed`: RepeatMasker (UCSC `rmsk.txt.gz`) elements of the transposon classes (`TE_CLASSES`), named
  `<repName>:<repFamily>:<repClass>`.
- `piRNA_cluster.bed`: piRNAclusterDB clusters from the UNITAS reference dump (`piRclDB`), located on this genome by
  unique exact matches of short anchors from both ends of each cluster sequence (Bowtie 1), so clusters published on
  another assembly (mouse: GRCm38) are placed on this one without a liftover; clusters that cannot be placed uniquely
  and whole are left out and counted. TEsmall matches strand only for miRNA and hairpin, so clusters get strand `.`.
- `rDNA.fa`: the NCBI rDNA unit (45S) given with `--rdna`, and the Ensembl rRNA and Mt_rRNA genes.
- `tDNA.fa`: GtRNAdb tRNA genes taken from this genome, sense strand, with the headers TEsmall's 3' tRF step parses,
  `<name>:<chrom>:<start>-<end>:<strand>` (BED start); each must equal GtRNAdb's own sequence.

Features on sequences the genome does not have (e.g. alternative contigs) are dropped and counted. A source that does
not fit the genome (a tRNA sequence that differs, a miRBase chromosome missing from the map) stops the build.

Usage: tesmall_db.py --genome FASTA --genome-index PREFIX --gtf GTF[.gz] --mirna-gff GFF3 --chr-map TSV
                     --trna-archive TAR.GZ --rmsk TXT.GZ --rdna FASTA --pirna-clusters FASTA
                     --annotation DIR --sequence DIR --record TSV [--source KEY=VALUE ...]
"""
import argparse
import collections
import gzip
import hashlib
import io
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mirflowz_chr_map import candidates  # noqa: E402

STRUCTURAL_BIOTYPES = ('rRNA', 'Mt_rRNA', 'Mt_tRNA', 'snRNA', 'snoRNA', 'scaRNA', 'misc_RNA', 'vault_RNA', 'ribozyme',
                       'sRNA', 'scRNA')
RDNA_BIOTYPES = ('rRNA', 'Mt_rRNA')
TE_CLASSES = ('LINE', 'SINE', 'LTR', 'DNA', 'RC', 'Retroposon')
CLASSES = ('structural_RNA', 'miRNA', 'hairpin', 'exon', 'TE', 'intron', 'piRNA_cluster')
ANCHOR = 30                       # anchor length (nt) at each end of a piRNA cluster
ANCHOR_OFFSETS = (0, 50, 100, 150, 200, 300, 500)
SPAN_TOLERANCE = 0.02             # placed cluster length may differ this much from its sequence (other assembly)


def fail(message):
    sys.exit('tesmall_db.py: ' + message)


def open_text(path):
    with open(path, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gzipped else open(path)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def genome_names(fai):
    """Sequence names and lengths from the genome's .fai."""
    with open(fai) as fh:
        return {line.split('\t')[0]: int(line.split('\t')[1]) for line in fh if line.strip()}


def rename(name, names):
    """The genome's name for a UCSC-style chromosome name (the mirflowz chromosome map rules), or None."""
    return next((c for c in candidates(name) if c in names), None)


def read_chr_map(path):
    """miRBase GFF chromosome -> genome name; names without target (removed contigs) map to None."""
    out = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if fields and fields[0]:
                out[fields[0]] = fields[1] if len(fields) > 1 and fields[1] else None
    return out


def attributes(text, sep='='):
    out = {}
    for part in text.strip().strip(';').split(';'):
        part = part.strip()
        if not part:
            continue
        if sep == '=':
            key, _, value = part.partition('=')
        else:
            key, _, value = part.partition(' ')
            value = value.strip().strip('"')
        out[key.strip()] = value
    return out


def mirna_beds(gff, chr_map):
    """(miRNA rows, hairpin rows) from a miRBase GFF3; a chromosome missing from the map stops."""
    mature, hairpin = [], []
    with open_text(gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] not in ('miRNA', 'miRNA_primary_transcript'):
                continue
            if f[0] not in chr_map:
                fail('{}: chromosome {} is not in the chromosome map'.format(gff, f[0]))
            chrom = chr_map[f[0]]
            if chrom is None:
                continue
            row = (chrom, int(f[3]) - 1, int(f[4]), attributes(f[8])['Name'], f[6])
            (mature if f[2] == 'miRNA' else hairpin).append(row)
    return mature, hairpin


def gtf_beds(gtf, names):
    """(structural rows, exon rows, intron rows, rDNA gene rows, dropped) from an Ensembl GTF."""
    structural, rdna, exons, dropped = [], [], collections.defaultdict(list), 0
    with open_text(gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('gene', 'exon'):
                continue
            if f[0] not in names:
                dropped += 1
                continue
            a = attributes(f[8], sep=' ')
            biotype = a.get('gene_biotype', '')
            gene = a.get('gene_name') or a['gene_id']
            start, end = int(f[3]) - 1, int(f[4])
            if f[2] == 'gene':
                if biotype in STRUCTURAL_BIOTYPES:
                    structural.append((f[0], start, end, '{}:{}'.format(biotype, gene), f[6]))
                if biotype in RDNA_BIOTYPES:
                    rdna.append((f[0], start, end, '{}:{}'.format(biotype, gene), f[6]))
            elif biotype not in STRUCTURAL_BIOTYPES and biotype != 'miRNA':
                exons[(a['transcript_id'], gene, f[0], f[6])].append((start, end))
    exon_rows, intron_rows = set(), set()
    for (_, gene, chrom, strand), spans in exons.items():
        spans.sort()
        exon_rows.update((chrom, s, e, gene, strand) for s, e in spans)
        intron_rows.update((chrom, e1, s2, gene, strand) for (_, e1), (s2, _) in zip(spans, spans[1:]) if s2 > e1)
    return structural, sorted(exon_rows), sorted(intron_rows), rdna, dropped


def rmsk_rows(path, names):
    """(TE rows, dropped) from UCSC rmsk.txt.gz."""
    rows, dropped = [], 0
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[11].rstrip('?') not in TE_CLASSES:
                continue
            chrom = rename(f[5], names)
            if chrom is None:
                dropped += 1
                continue
            rows.append((chrom, int(f[6]), int(f[7]), '{}:{}:{}'.format(f[10], f[12], f[11]), f[9]))
    return rows, dropped


def trna_records(archive):
    """GtRNAdb's high-confidence tRNA genes, `<db>-tRNAs.fa` of its tRNAs.tar.gz (the BED also lists low-scoring
    predictions and tRX loci): [(UCSC chrom, start0, end, name, strand, sequence)], name without the species prefix."""
    with tarfile.open(archive) as tar:
        members = {Path(m.name).name: m for m in tar.getmembers()}
        fa = [n for n in members if n.endswith('-tRNAs.fa') and '-mature-' not in n and '-filtered-' not in n]
        if len(fa) != 1:
            fail('{}: expected one <db>-tRNAs.fa, found {}'.format(archive, ', '.join(sorted(members))))
        fa_text = tar.extractfile(members[fa[0]]).read().decode()
    records = []
    for line in io.StringIO(fa_text):
        if line.startswith('>'):
            m = re.match(r'>(?:[A-Z][a-z]+_[a-z]+_)?(\S+) .* (\S+):(\d+)-(\d+) \(([+-])\)\s*$', line)
            if not m:
                fail('{}: unexpected header {}'.format(fa[0], line.strip()))
            records.append([m.group(2), int(m.group(3)) - 1, int(m.group(4)), m.group(1), m.group(5), []])
        elif records:
            records[-1][5].append(line.strip())
    return [tuple(r[:5]) + (''.join(r[5]).upper(),) for r in records]


def fetch(genome, regions):
    """Sequences of 1-based regions [(chrom, start, end, strand)] from an indexed FASTA (samtools faidx)."""
    out = {}
    for strand in ('+', '-'):
        todo = [r for r in regions if r[3] == strand]
        if not todo:
            continue
        with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as fh:
            fh.write(''.join('{}:{}-{}\n'.format(c, s, e) for c, s, e, _ in todo))
        try:
            cmd = ['samtools', 'faidx', str(genome), '-r', fh.name] + (['-i', '--mark-strand', 'no'] if strand == '-' else [])
            text = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        finally:
            os.unlink(fh.name)
        seqs = [''.join(block.splitlines()[1:]).upper() for block in text.split('>')[1:]]
        if len(seqs) != len(todo):
            fail('samtools faidx returned {} sequences for {} regions'.format(len(seqs), len(todo)))
        out.update(zip(todo, seqs))
    return out


def read_fasta(path):
    records, name = [], None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].strip()
                records.append([name, []])
            elif records:
                records[-1][1].append(line.strip())
    return [(n, ''.join(s)) for n, s in records]


def cluster_anchors(clusters):
    """[(anchor id, sequence)] for each cluster: anchors at ANCHOR_OFFSETS from both ends."""
    out = []
    for cid, seq in clusters:
        seq = seq.upper()
        for off in ANCHOR_OFFSETS:
            if len(seq) < 2 * (off + ANCHOR):
                break
            out.append(('{}|s|{}'.format(cid, off), seq[off:off + ANCHOR]))
            out.append(('{}|e|{}'.format(cid, off), seq[len(seq) - off - ANCHOR:len(seq) - off]))
    return out


def place_clusters(clusters, hits):
    """Place clusters from unique anchor hits {anchor id: (chrom, pos0, strand)}; returns (rows, reasons).

    The first start anchor and first end anchor (smallest offsets) that hit uniquely, on the same chromosome and in
    the same orientation, give the cluster's ends; the placed length must match the sequence within SPAN_TOLERANCE."""
    rows, reasons = [], collections.Counter()
    for cid, seq in clusters:
        n = len(seq)
        start = next(((off, hits[k]) for off in ANCHOR_OFFSETS if (k := '{}|s|{}'.format(cid, off)) in hits), None)
        end = next(((off, hits[k]) for off in ANCHOR_OFFSETS if (k := '{}|e|{}'.format(cid, off)) in hits), None)
        if start is None or end is None:
            reasons['no unique anchor at the {} of the cluster'.format('start' if start is None else 'end')] += 1
            continue
        (s_off, (s_chrom, s_pos, s_strand)), (e_off, (e_chrom, e_pos, e_strand)) = start, end
        if s_chrom != e_chrom or s_strand != e_strand:
            reasons['anchors on different chromosomes or strands'] += 1
            continue
        if s_strand == '+':
            begin, stop = s_pos - s_off, e_pos + ANCHOR + e_off
        else:
            begin, stop = e_pos - e_off, s_pos + ANCHOR + s_off
        if stop <= begin or abs((stop - begin) - n) > SPAN_TOLERANCE * n:
            reasons['placed length differs from the cluster sequence'] += 1
            continue
        rows.append((s_chrom, begin, stop, cid, '.'))
    return rows, reasons


def bowtie_hits(index, anchors):
    """Unique exact hits {anchor id: (chrom, pos0, strand)} of anchors on a Bowtie 1 index."""
    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / 'anchors.fa'
        fa.write_text(''.join('>{}\n{}\n'.format(a, s) for a, s in anchors))
        out = subprocess.run(['bowtie', '-f', '-v', '0', '-m', '1', '--best', '--strata', str(index), str(fa)],
                             check=True, capture_output=True, text=True).stdout
    hits = {}
    for line in out.splitlines():
        f = line.split('\t')
        hits[f[0]] = (f[2], int(f[3]), f[1])
    return hits


def pirna_clusters(path):
    """[(cluster id, sequence)] and {cluster id: declared (chrom, start0, end)} from a UNITAS piRclDB file."""
    clusters, declared = [], {}
    for header, seq in read_fasta(path):
        m = re.match(r'([A-Za-z]+_\d+)_(\S+):(\d+)\.\.(\d+)$', header.split()[0])
        if not m:
            fail('{}: unexpected header {}'.format(path, header))
        cid = 'piRNAclusterDB:' + m.group(1)
        clusters.append((cid, seq))
        declared[cid] = (m.group(2), int(m.group(3)) - 1, int(m.group(4)))
    return clusters, declared


def write_bed(path, rows):
    rows = sorted(set(rows))
    with open(path, 'w') as fh:
        fh.writelines('{}\t{}\t{}\t{}\t0\t{}\n'.format(c, s, e, n, st) for c, s, e, n, st in rows)
    return len(rows)


def write_fasta(path, records):
    with open(path, 'w') as fh:
        for name, seq in records:
            fh.write('>{}\n'.format(name))
            fh.writelines(seq[i:i + 60] + '\n' for i in range(0, len(seq), 60))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for opt in ('--genome', '--genome-index', '--gtf', '--mirna-gff', '--chr-map', '--trna-archive', '--rmsk', '--rdna',
                '--pirna-clusters', '--annotation', '--sequence', '--record'):
        parser.add_argument(opt, required=True, type=Path)
    parser.add_argument('--source', action='append', default=[], metavar='KEY=VALUE', help='extra record lines')
    args = parser.parse_args()
    names = genome_names(Path(str(args.genome) + '.fai'))
    args.annotation.mkdir(parents=True, exist_ok=True)
    args.sequence.mkdir(parents=True, exist_ok=True)
    record = [('sources', 'genome', str(args.genome))]

    mature, hairpin = mirna_beds(args.mirna_gff, read_chr_map(args.chr_map))
    structural, exons, introns, rdna_genes, gtf_dropped = gtf_beds(args.gtf, names)

    trnas, trna_dropped = [], 0
    for chrom, start, end, name, strand, seq in trna_records(args.trna_archive):
        target = rename(chrom, names)
        if target is None:
            trna_dropped += 1
            continue
        trnas.append((target, start, end, name, strand, seq))
    genomic = fetch(args.genome, [(c, s + 1, e, st) for c, s, e, _, st, _ in trnas])
    differ = [n for c, s, e, n, st, seq in trnas if genomic[(c, s + 1, e, st)] != seq]
    if differ:
        fail('{} GtRNAdb tRNA sequences differ from the genome (wrong assembly?), e.g. {}'.format(len(differ), ', '.join(differ[:5])))
    structural += [(c, s, e, 'tRNA:' + n, st) for c, s, e, n, st, _ in trnas]
    write_fasta(args.sequence / 'tDNA.fa', [('{}:{}:{}-{}:{}'.format(n, c, s, e, st), seq) for c, s, e, n, st, seq in trnas])

    rdna = read_fasta(args.rdna)
    genes = fetch(args.genome, [(c, s + 1, e, st) for c, s, e, _, st in rdna_genes])
    rdna += [('{}:{}:{}-{}'.format(n, c, s, e), genes[(c, s + 1, e, st)]) for c, s, e, n, st in rdna_genes]
    write_fasta(args.sequence / 'rDNA.fa', [(n.split()[0], s.upper()) for n, s in rdna])

    te, te_dropped = rmsk_rows(args.rmsk, names)

    clusters, declared = pirna_clusters(args.pirna_clusters)
    placed, reasons = place_clusters(clusters, bowtie_hits(args.genome_index, cluster_anchors(clusters)))
    agree = sum(1 for c, s, e, cid, _ in placed if (rename(declared[cid][0], names), declared[cid][1], declared[cid][2]) == (c, s, e))

    counts = {'miRNA': write_bed(args.annotation / 'miRNA.bed', mature),
              'hairpin': write_bed(args.annotation / 'hairpin.bed', hairpin),
              'structural_RNA': write_bed(args.annotation / 'structural_RNA.bed', structural),
              'exon': write_bed(args.annotation / 'exon.bed', exons),
              'intron': write_bed(args.annotation / 'intron.bed', introns),
              'TE': write_bed(args.annotation / 'TE.bed', te),
              'piRNA_cluster': write_bed(args.annotation / 'piRNA_cluster.bed', placed)}
    for path in (args.gtf, args.mirna_gff, args.chr_map, args.trna_archive, args.rmsk, args.rdna, args.pirna_clusters):
        record.append(('sha256', path.name, sha256(path)))
    record += [('features', cls, counts[cls]) for cls in CLASSES]
    record += [('sequences', 'tDNA', len(trnas)), ('sequences', 'rDNA', len(rdna)),
               ('dropped', 'GTF features on sequences not in the genome', gtf_dropped),
               ('dropped', 'GtRNAdb tRNAs on sequences not in the genome', trna_dropped),
               ('dropped', 'RepeatMasker elements on sequences not in the genome', te_dropped),
               ('piRNA clusters', 'in piRclDB', len(clusters)), ('piRNA clusters', 'placed', len(placed)),
               ('piRNA clusters', 'placed at their published coordinates', agree)]
    record += [('piRNA clusters', 'not placed: ' + reason, n) for reason, n in sorted(reasons.items())]
    record += [('source', *item.split('=', 1)) for item in args.source]
    with open(args.record, 'w') as fh:
        fh.write('section\tkey\tvalue\n')
        fh.writelines('{}\t{}\t{}\n'.format(*row) for row in record)


if __name__ == '__main__':
    main()
