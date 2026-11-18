"""Small-RNA filtering with cutadapt only (filter.trim.quantifier: cutadapt, #350).

Three cutadapt passes per sample:

1. the 3' adapter (`adapter`) with a minimum overlap of THREE_PRIME_MIN_OVERLAP (cutadapt aligns with indels);
2. for kits with random ends (`random_ends`, e.g. 4 for NEXTflex v3) that many bases cut from both ends, after the
   adapter search (cutadapt applies -u before adapter search, and cutting first truncates adapter dimers whose adapter
   starts within the first bases); reads containing the 5' adapter discarded; the kit's minimum length and N limit;
3. reads still containing 3' adapter or Illumina primer sequence discarded (#353): mostly inserts whose 3' adapter has
   too many sequencing errors for pass 1 (in public plasma data often the RNY4 Y RNA fragment), adapter dimers with an
   indel or duplication in the adapter, and reads of primer sequence only. With that much adapter left, no quantifier
   annotates them (UNITAS: "no annotation"), so discarding them changes no count; it names them in the report.
   Reads of the kit's reagent oligos (`reagent_oligos`, #380) are discarded in the same pass and counted apart: these
   do change counts (UNITAS counted the TruSeq Stop Oligo as protein-coding reads).

Pass 3 is a separate pass so that its `discard_trimmed` counts exactly these reads: cutadapt matches adapters before it
applies the length and N filters, so per-record match counts of a combined pass include reads those filters removed.
Its records (three_prime_artefact_fasta) are every piece of at least THREE_PRIME_ARTEFACT_PIECE nt of the kit's 3'
adapter followed by the start of Illumina's index-read primer (ILLUMINA_INDEX_PRIMER_START, as in the read:
`...TGGAATTCTCGGGTGCCAAGG` `AACTCCAGTCAC` `<index>` `ATCTCGTATGCCGTCTTCTGCTTG`), and of Illumina's P7 sequence
(ILLUMINA_P7); the index between them differs per sample and is not used. The Illumina sequences are the same for every
kit, so they are constants rather than kit settings. A read with both a 5' artefact and 3' artefact sequence counts as
a 5' artefact (pass 2 removes it first). The 16-nt minimum is essential: with 14-nt pieces (10% errors), pieces at the
adapter's end matched human 18S rRNA (around position 400), hsa-miR-2116-3p and a misc_RNA; from 16 nt no human mature
miRNA (miRBase 23.1) or Ensembl/GtRNAdb/SILVA tRNA, rRNA, snRNA, snoRNA, scaRNA, vault or misc_RNA sequence matches.

The kit's reagent oligos follow (reagent_oligo_fasta), two records each:
- the oligo as a 5' adapter with a minimum overlap of REAGENT_OLIGO_MIN_OVERLAP nt and FIVE_PRIME_MAX_ERROR_RATE: the
  oligo complete anywhere in the read, or its 3' end of at least that length at the read start, as in the 5'-truncated
  copies of the TruSeq Small RNA Stop Oligo (STP), which pass the other filters and which UNITAS counts as
  protein-coding reads (COL5A3). Its one allowed error keeps copies with a sequencing error; the cost is that a
  fragment of one of those 78 transcripts that starts exactly where it matches the oligo's end is discarded too;
- its last REAGENT_OLIGO_MIN_OVERLAP nt, exactly, anywhere in the read: copies after a few other bases, with a deletion
  near their start or repeated end to end (0.27% of the reads the first record kept in a public TruSeq plasma
  library, e.g. AA + CACCACGTTCCCGTGG). Exactly, because with one error the STP's end matches 78 human transcripts
  (Ensembl 116 cDNA, COL5A3 among them); exactly, it matches one, and no human mature miRNA either way.
Chimeras of an oligo and an insert are discarded too (7 of 5.8 million reads in public TruSeq plasma data at 14 nt).
Without the setting (or with null) a kit has none, as with random_ends. A read that matches both a 3' artefact and an
oligo record counts for the record cutadapt finds best; merge_cutadapt.py counts the two kinds apart by record name.

5' adapter artefacts are fragments of the 5' adapter (`adapter5`) that start anywhere in it, usually end before its
end, and sometimes follow a few other bases. cutadapt accepts a partial 5' adapter only when it runs to the adapter's 3'
end, and finds an internal 5' adapter only when it is complete, so two sets are given (five_prime_fasta):

- the adapter truncated at every end position (adapter5[:j]) as non-internal adapters (`X` prefix), minimum overlap
  FIVE_PRIME_MIN_OVERLAP: a fragment adapter5[i:j] at the read start is a partial match of adapter5[:j], and these
  records never match inside the read;
- every piece of the adapter of at least `five_prime_piece` nt (a kit setting; adapter5[i:i+L], prefixes included),
  which matches when complete anywhere in the read, e.g. after junk bases (`AAACAGTCCGACGATC`).

Every record has the error rate FIVE_PRIME_MAX_ERROR_RATE set explicitly (cutadapt's default, mismatches and indels):
one error in matches of 12-19 nt, two in matches of 20-26 nt. A count of one error is not possible: cutadapt scales an
absolute count to the full record, so partial matches would get none; with one error, or without indels, real artefacts
(e.g. adapter fragments with a two-base deletion) were kept. The lengths are essential: a 10-nt overlap also removed
miR-22-5p and miR-139-5p (their first 11 bases match the adapter with one mismatch), and 12-13 nt pieces anywhere match
miR-625-5p; with 12 at the read start and 14 anywhere, the rule matches only miR-3168-5p among the 2,789 human mature
miRNAs of miRBase 23.1. The piece length is a kit setting because it depends on the adapter: 14 for the TruSeq/NEXTflex
v3 adapter, 15 for the NEXTflex v4 adapter (TCTTTCCCTACACGACGCTCTTCCGATCT), whose 14-nt pieces match hsa-miR-6865-3p
with one error; from 15 nt nothing in human mature miRNAs or the small-RNA references matches. 15 nt for the v3/TruSeq
adapter would keep artefacts (in a low-yield public plasma library 0.3% of kept reads, mostly 11-13 nt of the adapter's
end after a few other bases, e.g. AAAAGTCCGACGATC) (#350). The match only flags the read (--action=none), so the report counts discarded artefacts as such
rather than as too short.

Included by smallrna/rules/filter/cutadapt.smk; also a command-line script for the adapter FASTAs:
    cutadapt_config.py five-prime-fasta ADAPTER5 PIECE > five_prime_adapters.fa
    cutadapt_config.py artefact-fasta ADAPTER3 [OLIGO ...] > artefacts.fa
"""
import argparse
import re

THREE_PRIME_MIN_OVERLAP = 5
FIVE_PRIME_MIN_OVERLAP = 12
MIN_FIVE_PRIME_PIECE = 14  # 12-13-nt pieces anywhere match miR-625-5p
FIVE_PRIME_MAX_ERROR_RATE = 0.1
THREE_PRIME_ARTEFACT_PIECE = 16
THREE_PRIME_ARTEFACT_RECORD, REAGENT_OLIGO_RECORD = 'three_prime_artefact_', 'reagent_oligo_'  # record name prefixes, read by merge_cutadapt.py
REAGENT_OLIGO_MIN_OVERLAP = 14  # in TruSeq plasma: 7 STP-miRNA chimeras of 5.8 M reads discarded at 14 nt, 43 at 10 (#380)
ILLUMINA_INDEX_PRIMER_START = 'AACTCCAGTCAC'
ILLUMINA_P7 = 'ATCTCGTATGCCGTCTTCTGCTTG'
DNA = re.compile(r'^[ACGT]+\Z')  # \Z, not $: '$' also matches before a trailing newline


def check_oligos(oligos):
    """Fail unless oligos is a list of distinct A/C/G/T sequences of at least REAGENT_OLIGO_MIN_OVERLAP nt."""
    if not isinstance(oligos, list) or not all(isinstance(o, str) and DNA.match(o) and len(o) >= REAGENT_OLIGO_MIN_OVERLAP for o in oligos):
        raise ValueError('a list of A/C/G/T sequences of at least {} nt, got {!r}'.format(REAGENT_OLIGO_MIN_OVERLAP, oligos))
    ends = [o[-REAGENT_OLIGO_MIN_OVERLAP:] for o in oligos]
    if len(set(ends)) != len(ends):
        raise ValueError('oligos must differ in their last {} nt (an oligo listed twice?): {!r}'.format(REAGENT_OLIGO_MIN_OVERLAP, oligos))


def cutadapt_settings(config):
    """Validated settings from config['adapter'] and config['filter']['trim']['cutadapt']; fails on anything missing
    except the optional random_ends (0) and reagent_oligos (none)."""
    settings = (config.get('filter', {}).get('trim', {}) or {}).get('cutadapt') or {}
    kit = config.get('libprepkit', 'this kit')
    if not settings:
        raise ValueError('filter.trim.quantifier: cutadapt needs filter.trim.cutadapt settings, none for {!r}; '
                         'set them in libprep.config or the project configuration'.format(kit))
    unknown = sorted(set(settings) - {'adapter5', 'five_prime_piece', 'random_ends', 'minimum_length', 'max_n', 'max_ee', 'reagent_oligos'})
    if unknown:
        raise ValueError('filter.trim.cutadapt: unknown setting(s) {}'.format(', '.join(unknown)))
    adapter3 = config.get('adapter')
    if not isinstance(adapter3, str) or not DNA.match(adapter3):
        raise ValueError('filter.trim.quantifier: cutadapt needs the kit\'s 3\' adapter as an A/C/G/T sequence in '
                         '`adapter`, got {!r}'.format(adapter3))
    adapter5 = settings.get('adapter5')
    if not isinstance(adapter5, str) or not DNA.match(adapter5) or len(adapter5) < FIVE_PRIME_MIN_OVERLAP:
        raise ValueError('filter.trim.cutadapt.adapter5: the kit\'s 5\' adapter as an A/C/G/T sequence of at least {} nt, '
                         'got {!r}'.format(FIVE_PRIME_MIN_OVERLAP, adapter5))
    piece = settings.get('five_prime_piece')
    if not isinstance(piece, int) or isinstance(piece, bool) or not MIN_FIVE_PRIME_PIECE <= piece <= len(adapter5):
        raise ValueError('filter.trim.cutadapt.five_prime_piece: the minimum length of 5\' adapter pieces matched anywhere, '
                         'from {} nt to the adapter length, got {!r}'.format(MIN_FIVE_PRIME_PIECE, piece))
    random_ends = settings.get('random_ends', 0)
    minimum_length = settings.get('minimum_length')
    max_n = settings.get('max_n')
    max_ee = settings.get('max_ee')
    if not isinstance(random_ends, int) or isinstance(random_ends, bool) or random_ends < 0:
        raise ValueError('filter.trim.cutadapt.random_ends: a number of bases >= 0, got {!r}'.format(random_ends))
    if not isinstance(minimum_length, int) or isinstance(minimum_length, bool) or minimum_length < 1:
        raise ValueError('filter.trim.cutadapt.minimum_length: the final minimum read length (>= 1), got {!r}'.format(minimum_length))
    if max_n is not None and (not isinstance(max_n, int) or isinstance(max_n, bool) or max_n < 0):
        raise ValueError('filter.trim.cutadapt.max_n: a number of N bases >= 0, or null, got {!r}'.format(max_n))
    if max_ee is not None and (isinstance(max_ee, bool) or not isinstance(max_ee, (int, float)) or max_ee <= 0):
        raise ValueError('filter.trim.cutadapt.max_ee: expected errors > 0, or null, got {!r}'.format(max_ee))
    oligos = settings.get('reagent_oligos')
    oligos = [] if oligos is None else oligos
    try:
        check_oligos(oligos)
    except ValueError as error:
        raise ValueError('filter.trim.cutadapt.reagent_oligos: {}'.format(error))
    return {'adapter3': adapter3, 'adapter5': adapter5, 'five_prime_piece': piece, 'random_ends': random_ends, 'minimum_length': minimum_length,
            'max_n': max_n, 'max_ee': max_ee, 'reagent_oligos': oligos}


def adapter_args(settings):
    """Options of the first pass: the 3' adapter."""
    return "-a '{};min_overlap={}'".format(settings['adapter3'], THREE_PRIME_MIN_OVERLAP)


def filter_args(settings):
    """Options of the second pass, except the 5' adapter FASTA (-g file:...)."""
    args = []
    if settings['random_ends']:
        args += ['-u {0} -u -{0}'.format(settings['random_ends'])]
    args += ['--action=none --discard-trimmed', '-m {}'.format(settings['minimum_length'])]
    if settings['max_n'] is not None:
        args += ['--max-n {}'.format(settings['max_n'])]
    if settings['max_ee'] is not None:
        args += ['--max-ee {:g}'.format(settings['max_ee'])]
    return ' '.join(args)


def five_prime_fasta(adapter5, piece):
    """cutadapt adapter FASTA: adapter5 truncated at every end position (from FIVE_PRIME_MIN_OVERLAP nt) as non-internal
    adapters, then every piece of at least `piece` nt; all with FIVE_PRIME_MAX_ERROR_RATE."""
    truncated = [adapter5[:j] for j in range(FIVE_PRIME_MIN_OVERLAP, len(adapter5) + 1)]
    pieces = sorted({adapter5[i:i + n] for n in range(piece, len(adapter5) + 1) for i in range(len(adapter5) - n + 1)},
                    key=lambda p: (-len(p), adapter5.index(p)))
    rate = FIVE_PRIME_MAX_ERROR_RATE
    return (''.join('>five_prime_{0}\nX{1};min_overlap={2};max_error_rate={3}\n'.format(len(t), t, FIVE_PRIME_MIN_OVERLAP, rate)
                    for t in truncated) +
            ''.join('>five_prime_piece_{0}_{1}\n{2};min_overlap={3};max_error_rate={4}\n'.format(
                adapter5.index(p), len(p), p, piece, rate) for p in pieces))


def artefact_args():
    """Options of the third pass, except the artefact FASTA (-g file:...): only flag and discard. cutadapt's default of
    one match per read (-n 1) lets merge_cutadapt.py split the discarded reads by record; keep it."""
    return '--action=none --discard-trimmed'


def three_prime_artefact_fasta(adapter3):
    """cutadapt adapter FASTA of the third pass: every piece of at least THREE_PRIME_ARTEFACT_PIECE nt of adapter3 plus
    ILLUMINA_INDEX_PRIMER_START, and of ILLUMINA_P7, matched when complete anywhere in the read, with
    FIVE_PRIME_MAX_ERROR_RATE."""
    records = []
    for tag, sequence in (('adapter', adapter3 + ILLUMINA_INDEX_PRIMER_START), ('p7', ILLUMINA_P7)):
        pieces = sorted({sequence[i:i + n] for n in range(THREE_PRIME_ARTEFACT_PIECE, len(sequence) + 1)
                         for i in range(len(sequence) - n + 1)}, key=lambda p: (-len(p), sequence.index(p)))
        records += ['>{}{}_{}_{}\n{};min_overlap={};max_error_rate={}\n'.format(
            THREE_PRIME_ARTEFACT_RECORD, tag, sequence.index(p), len(p), p, THREE_PRIME_ARTEFACT_PIECE, FIVE_PRIME_MAX_ERROR_RATE) for p in pieces]
    return ''.join(records)


def reagent_oligo_fasta(oligos):
    """cutadapt adapter FASTA records of the kit's reagent oligos: the oligo as a 5' adapter with
    REAGENT_OLIGO_MIN_OVERLAP and FIVE_PRIME_MAX_ERROR_RATE (complete anywhere, or its 3' end at the read start), and
    its last REAGENT_OLIGO_MIN_OVERLAP nt without errors, anywhere (left out for an oligo of just that length)."""
    n, records = REAGENT_OLIGO_MIN_OVERLAP, []
    for i, o in enumerate(oligos, 1):
        records.append('>{}{}_{}\n{};min_overlap={};max_error_rate={}\n'.format(REAGENT_OLIGO_RECORD, i, len(o), o, n, FIVE_PRIME_MAX_ERROR_RATE))
        if len(o) > n:
            records.append('>{}{}_end_{}\n{};min_overlap={};max_error_rate=0\n'.format(REAGENT_OLIGO_RECORD, i, n, o[-n:], n))
    return ''.join(records)


def artefact_fasta(adapter3, oligos):
    """cutadapt adapter FASTA of the third pass: the 3' artefact records, then the reagent oligos."""
    return three_prime_artefact_fasta(adapter3) + reagent_oligo_fasta(oligos)


def main():
    parser = argparse.ArgumentParser(description='Write the adapter FASTAs for cutadapt.')
    sub = parser.add_subparsers(dest='command', required=True)
    fasta = sub.add_parser('five-prime-fasta')
    fasta.add_argument('adapter5')
    fasta.add_argument('piece', type=int)
    artefacts = sub.add_parser('artefact-fasta')
    artefacts.add_argument('adapter3')
    artefacts.add_argument('oligos', nargs='*', help="the kit's reagent oligos")
    args = parser.parse_args()
    if args.command == 'artefact-fasta':
        if not DNA.match(args.adapter3):
            parser.error('adapter3 must be an A/C/G/T sequence')
        try:
            check_oligos(args.oligos)
        except ValueError as error:
            parser.error('reagent oligos: {}'.format(error))
        print(artefact_fasta(args.adapter3, args.oligos), end='')
        return
    if not DNA.match(args.adapter5) or len(args.adapter5) < FIVE_PRIME_MIN_OVERLAP:
        parser.error('adapter5 must be an A/C/G/T sequence of at least {} nt'.format(FIVE_PRIME_MIN_OVERLAP))
    if not MIN_FIVE_PRIME_PIECE <= args.piece <= len(args.adapter5):
        parser.error('piece must be from {} nt to the adapter length'.format(MIN_FIVE_PRIME_PIECE))
    print(five_prime_fasta(args.adapter5, args.piece), end='')


if __name__ == '__main__':
    main()
