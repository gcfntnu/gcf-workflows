#-*- mode:snakemake -*-
"""Spike-in sets of `filter.spikein.ref`, shared by the bowtie, UNITAS and miRge3 routes.

Spike-ins are a library prep choice. `filter.spikein.ref` lists sets, comma-separated: gcfdb sets by name (sequences
in `EXT_DIR/spikein/fasta/<set>.fa`), or a project's own FASTA by its path (ending in .fa or .fasta, relative to the
working directory), e.g. synthetic non-human sequences added to the libraries. Whatever the route, the same sequences
are used: the bowtie route counts each spike-in (rules/filter/spikein_bowtie.smk), UNITAS adds them to its reference
with `-refseq`, and miRge3 gets a spike-in index built from them.

`filter.spikein.quantifier`: `skip`, `bowtie` or `unitas`. `filter.spikein.filter: true` removes the reads matching a
spike-in before quantification; bowtie route only.
"""
import re

SPIKEIN_SETS = ('ERCC92', 'smallrna_calibrators')
SPIKEIN_QUANTIFIERS = ('skip', 'bowtie', 'unitas')
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

if config['filter']['spikein'].get('quantifier') not in SPIKEIN_QUANTIFIERS:
    raise ValueError('filter.spikein.quantifier: smallrna supports {}, got {!r}'.format(
        ', '.join(SPIKEIN_QUANTIFIERS), config['filter']['spikein'].get('quantifier')))
if config['filter']['spikein']['quantifier'] != 'skip' and not SPIKEIN_REF_SETS:
    raise ValueError('filter.spikein.quantifier is {} but filter.spikein.ref lists no spike-in sets'.format(
        config['filter']['spikein']['quantifier']))
if config['filter']['spikein'].get('filter') not in (True, False):
    raise ValueError('filter.spikein.filter: true or false, got {!r}'.format(config['filter']['spikein'].get('filter')))
if config['filter']['spikein']['filter'] and config['filter']['spikein']['quantifier'] != 'bowtie':
    raise ValueError('filter.spikein.filter: true removes spike-in reads with the bowtie route only; set '
                     'filter.spikein.quantifier: bowtie, got {!r}'.format(config['filter']['spikein']['quantifier']))


def spikein_fasta():
    """FASTA files of the configured spike-in sets: gcfdb's for set names, the project's own for paths."""
    return [join(EXT_DIR, 'spikein', 'fasta', name + '.fa') if name in SPIKEIN_SETS else name for name in SPIKEIN_REF_SETS]
