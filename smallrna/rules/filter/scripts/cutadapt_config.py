"""Small-RNA filtering with cutadapt only (filter.trim.quantifier: cutadapt, #350).

Two cutadapt passes per sample:

1. the 3' adapter (`adapter`) with a minimum overlap of THREE_PRIME_MIN_OVERLAP (cutadapt aligns with indels);
2. for kits with random ends (`random_ends`, e.g. 4 for NEXTflex v3) that many bases cut from both ends, after the
   adapter search (cutadapt applies -u before adapter search, and cutting first truncates adapter dimers whose adapter
   starts within the first bases); reads containing the 5' adapter discarded; the kit's minimum length and N limit.

5' adapter artefacts are fragments of the 5' adapter (`adapter5`) that start anywhere in it, usually end before its
end, and sometimes follow a few other bases. cutadapt accepts a partial 5' adapter only when it runs to the adapter's 3'
end, and finds an internal 5' adapter only when it is complete, so two sets are given (five_prime_fasta):

- the adapter truncated at every end position (adapter5[:j]), minimum overlap FIVE_PRIME_MIN_OVERLAP: a fragment
  adapter5[i:j] at the read start is a partial match of adapter5[:j];
- every piece of the adapter of at least FIVE_PRIME_PIECE nt (adapter5[i:i+L]), which matches when complete anywhere in
  the read, e.g. after junk bases (`AAACAGTCCGACGATC`).

With cutadapt's default error rate both allow one mismatch. The lengths are essential: a 10-nt overlap also removed
miR-22-5p and miR-139-5p (their first 11 bases match the adapter with one mismatch), and 12-13 nt pieces anywhere match
miR-625-5p; with 12 at the read start and 14 anywhere, the rule matches only miR-3168-5p among the 2,789 human mature
miRNAs of miRBase 23.1. The match only flags the read (--action=none), so the report counts discarded artefacts as such
rather than as too short.

Included by smallrna/rules/filter/cutadapt.smk; also a command-line script for the truncated adapter FASTA:
    cutadapt_config.py five-prime-fasta ADAPTER5 > five_prime_adapters.fa
"""
import argparse
import re

THREE_PRIME_MIN_OVERLAP = 5
FIVE_PRIME_MIN_OVERLAP = 12
FIVE_PRIME_PIECE = 14
DNA = re.compile(r'^[ACGT]+$')


def cutadapt_settings(config):
    """Validated settings from config['adapter'] and config['filter']['trim']['cutadapt']; fails on anything missing."""
    settings = (config.get('filter', {}).get('trim', {}) or {}).get('cutadapt') or {}
    kit = config.get('libprepkit', 'this kit')
    if not settings:
        raise ValueError('filter.trim.quantifier: cutadapt needs filter.trim.cutadapt settings, none for {!r}; '
                         'set them in libprep.config or the project configuration'.format(kit))
    unknown = sorted(set(settings) - {'adapter5', 'random_ends', 'minimum_length', 'max_n', 'max_ee'})
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
    return {'adapter3': adapter3, 'adapter5': adapter5, 'random_ends': random_ends, 'minimum_length': minimum_length,
            'max_n': max_n, 'max_ee': max_ee}


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


def five_prime_fasta(adapter5):
    """cutadapt adapter FASTA: adapter5 truncated at every end position (from FIVE_PRIME_MIN_OVERLAP nt), then every
    other piece of at least FIVE_PRIME_PIECE nt."""
    truncated = [adapter5[:j] for j in range(FIVE_PRIME_MIN_OVERLAP, len(adapter5) + 1)]
    pieces = sorted({adapter5[i:i + n] for n in range(FIVE_PRIME_PIECE, len(adapter5) + 1)
                     for i in range(len(adapter5) - n + 1)} - set(truncated), key=lambda p: (-len(p), adapter5.index(p)))
    return (''.join('>five_prime_{0}\n{1};min_overlap={2}\n'.format(len(t), t, FIVE_PRIME_MIN_OVERLAP) for t in truncated) +
            ''.join('>five_prime_piece_{0}_{1}\n{2};min_overlap={3}\n'.format(adapter5.index(p), len(p), p, FIVE_PRIME_PIECE)
                    for p in pieces))


def main():
    parser = argparse.ArgumentParser(description='Write the truncated 5\' adapter FASTA for cutadapt.')
    sub = parser.add_subparsers(dest='command', required=True)
    fasta = sub.add_parser('five-prime-fasta')
    fasta.add_argument('adapter5')
    args = parser.parse_args()
    if not DNA.match(args.adapter5) or len(args.adapter5) < FIVE_PRIME_MIN_OVERLAP:
        parser.error('adapter5 must be an A/C/G/T sequence of at least {} nt'.format(FIVE_PRIME_MIN_OVERLAP))
    print(five_prime_fasta(args.adapter5), end='')


if __name__ == '__main__':
    main()
