#-*- mode:snakemake -*-
"""Spike-in sets of `filter.spikein.ref`, shared by the bowtie, UNITAS and miRge3 routes.

Spike-ins are set per project, in its own config, not per library prep kit, and `filter.spikein.ref` alone turns them
on: with no sets listed (the default, null) no spike-in step runs. It lists sets, comma-separated: gcfdb sets by name
(sequences in `EXT_DIR/spikein/fasta/<set>.fa`), or a project's own FASTA by its path (ending in .fa or .fasta,
relative to the working directory), e.g. synthetic non-human sequences added to the libraries.

`filter.spikein.quantifier`, optional, says who estimates them, from the same sequences:
- `bowtie` (default): a prefilter (rules/filter/spikein_bowtie.smk) counts each spike-in in the filtered reads and
  removes those reads, so every later step reads reads without spike-ins and none estimates them: UNITAS gets no
  spike-in reference and miRge3 no spike-in index;
- `native`: every quantifier that supports spike-ins estimates them itself from the reads with the spike-ins: UNITAS
  adds them to its reference with `-refseq`, miRge3 gets a spike-in index built from them. The bowtie prefilter still
  runs for the steps without spike-in support (mirflowz, miRTrace, STAR, FastQC), and its counts per spike-in are
  delivered as well. Setting it without sets stops the workflow.
Former settings stop the workflow: `skip` (no spike-ins is now an empty `filter.spikein.ref`), `unitas` (now `native`)
and `filter.spikein.filter` (now implied by `bowtie`).
A project's own FASTA is used as a copy in DNA letters (U to T), so RNA-letter files work for every route.
"""
import os
import re

SPIKEIN_SETS = ('ERCC92', 'smallrna_calibrators')
SPIKEIN_QUANTIFIERS = ('bowtie', 'native')
CUSTOM_SPIKEIN = re.compile(r'\.(fa|fasta)$')


def check_custom_spikein_fasta(path):
    """A project's spike-in FASTA: readable, at least one record, every record named, with an A/C/G/T/U sequence."""
    try:
        with open(path) as fh:
            lines = [line.strip() for line in fh if line.strip()]
    except OSError as error:
        raise ValueError('filter.spikein.ref: {!r} is neither a gcfdb spike-in set ({}) nor a readable FASTA: {}'.format(
            path, ', '.join(SPIKEIN_SETS), error))
    records = []
    for line in lines:
        if line.startswith('>'):
            if not line[1:].split():
                raise ValueError('filter.spikein.ref: {}: a record without a name'.format(path))
            records.append([line[1:].split()[0], ''])
        elif not records:
            raise ValueError('filter.spikein.ref: {}: not a FASTA file (sequence before the first header)'.format(path))
        elif not re.match(r'^[ACGTUacgtu]+$', line):
            raise ValueError('filter.spikein.ref: {}: record {!r} has bases other than A/C/G/T/U'.format(path, records[-1][0]))
        else:
            records[-1][1] += line
    if not records:
        raise ValueError('filter.spikein.ref: {}: no FASTA records'.format(path))
    empty = [name for name, sequence in records if not sequence]
    if empty:
        raise ValueError('filter.spikein.ref: {}: record {!r} has no sequence'.format(path, empty[0]))


def spikein_sets():
    """The comma-separated spike-in sets in filter.spikein.ref, in the given order; empty when none."""
    ref = config['filter']['spikein'].get('ref')
    names = [name.strip() for name in str(ref or '').split(',') if name.strip()]
    unknown = [name for name in names if name not in SPIKEIN_SETS and not CUSTOM_SPIKEIN.search(name)]
    if unknown:
        raise ValueError('filter.spikein.ref: unsupported value {!r}; use a comma-separated list of {}, or paths to the '
                         'project\'s own FASTA files (.fa, .fasta)'.format(ref, ', '.join(SPIKEIN_SETS)))
    if len(set(names)) != len(names):
        raise ValueError('filter.spikein.ref: a set is listed twice in {!r}'.format(ref))
    for name in names:
        if name not in SPIKEIN_SETS:
            check_custom_spikein_fasta(name)
    return names

SPIKEIN_REF_SETS = spikein_sets()

SPIKEIN_QUANTIFIER = config['filter']['spikein'].get('quantifier')
if SPIKEIN_QUANTIFIER == 'unitas':
    raise ValueError('filter.spikein.quantifier: unitas is now called native (every quantifier that supports spike-ins '
                     'estimates them: UNITAS with -refseq, miRge3 with a spike-in index); set filter.spikein.quantifier: native')
if SPIKEIN_QUANTIFIER == 'skip':
    raise ValueError('filter.spikein.quantifier: skip is no longer used: spike-ins are off when filter.spikein.ref lists '
                     'no sets; remove filter.spikein.quantifier')
if SPIKEIN_QUANTIFIER not in SPIKEIN_QUANTIFIERS:
    raise ValueError('filter.spikein.quantifier: smallrna supports {}, got {!r}'.format(', '.join(SPIKEIN_QUANTIFIERS), SPIKEIN_QUANTIFIER))
if SPIKEIN_QUANTIFIER == 'native' and not SPIKEIN_REF_SETS:
    raise ValueError('filter.spikein.quantifier is native but filter.spikein.ref lists no spike-in sets')
if 'filter' in config['filter']['spikein']:
    raise ValueError('filter.spikein.filter is no longer used: filter.spikein.quantifier: bowtie counts and removes spike-in '
                     'reads before quantification, native leaves them to the quantifiers; remove filter.spikein.filter')

# Spike-ins are on exactly when sets are listed; native additionally hands them to the quantifiers that estimate them.
SPIKEIN_ACTIVE = bool(SPIKEIN_REF_SETS)
SPIKEIN_NATIVE = SPIKEIN_ACTIVE and SPIKEIN_QUANTIFIER == 'native'

# Switching spike-ins off in a project that had them: the report reads all of bfq/figs and delivery takes bfq/exprs, so
# outputs of an earlier run would still be reported and delivered. Stop rather than delete delivery files silently.
SPIKEIN_BFQ_OUTPUTS = (join(BFQ_INTERIM, 'figs', 'spikein_mqc.yaml'), join(BFQ_INTERIM, 'exprs', 'spikein_counts.tsv'))
if not SPIKEIN_ACTIVE:
    SPIKEIN_STALE = [path for path in SPIKEIN_BFQ_OUTPUTS if os.path.lexists(path)]
    if SPIKEIN_STALE:
        raise ValueError('filter.spikein.ref lists no spike-in sets, but spike-in outputs from an earlier run remain and '
                         'would be reported and delivered: {}; delete them, or list the sets again'.format(', '.join(SPIKEIN_STALE)))


# A project's own FASTA files, by their name in the filter interim directory: rule spikein_project_fasta
# (filter/spikein_bowtie.smk) writes each in DNA letters there, so every route reads the same sequences.
SPIKEIN_CUSTOM = {'custom{}'.format(i): path for i, path in enumerate([n for n in SPIKEIN_REF_SETS if n not in SPIKEIN_SETS], 1)}


def native_spikein_fasta():
    """The spike-in FASTA files for quantifiers with their own spike-in support (UNITAS, miRge3): those of
    spikein_fasta() with filter.spikein.quantifier: native, none otherwise."""
    return spikein_fasta() if SPIKEIN_NATIVE else []


def spikein_fasta():
    """FASTA files of the configured spike-in sets, in order: gcfdb's for set names, the DNA copy for project files."""
    custom = {path: name for name, path in SPIKEIN_CUSTOM.items()}
    return [join(EXT_DIR, 'spikein', 'fasta', name + '.fa') if name in SPIKEIN_SETS
            else join(FILTER_INTERIM, 'spikein', 'sets', custom[name] + '.fa') for name in SPIKEIN_REF_SETS]
