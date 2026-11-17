#-*- mode:snakemake -*-
"""Spike-in sets of `filter.spikein.ref`, shared by the bowtie, UNITAS and miRge3 routes.

Spike-ins are a library prep choice. `filter.spikein.ref` lists sets, comma-separated: gcfdb sets by name (sequences
in `EXT_DIR/spikein/fasta/<set>.fa`), or a project's own FASTA by its path (ending in .fa or .fasta, relative to the
working directory), e.g. synthetic non-human sequences added to the libraries.

`filter.spikein.quantifier` says who estimates them, from the same sequences:
- `bowtie`: a prefilter (rules/filter/spikein_bowtie.smk) counts each spike-in in the filtered reads and removes those
  reads, so every later step reads reads without spike-ins and none estimates them: UNITAS gets no spike-in reference
  and miRge3 no spike-in index;
- `native`: every quantifier that supports spike-ins estimates them itself from the reads with the spike-ins: UNITAS
  adds them to its reference with `-refseq`, miRge3 gets a spike-in index built from them. The bowtie prefilter still
  runs for the steps without spike-in support (mirflowz, miRTrace, STAR, FastQC), and its counts per spike-in are
  delivered as well;
- `skip`: no spike-in estimation; then no sets may be listed.
`unitas`, the former name of `native`, and `filter.spikein.filter`, now implied by `bowtie`, stop the workflow.
A project's own FASTA is used as a copy in DNA letters (U to T), so RNA-letter files work for every route.
"""
import re

SPIKEIN_SETS = ('ERCC92', 'smallrna_calibrators')
SPIKEIN_QUANTIFIERS = ('skip', 'bowtie', 'native')
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
if SPIKEIN_QUANTIFIER not in SPIKEIN_QUANTIFIERS:
    raise ValueError('filter.spikein.quantifier: smallrna supports {}, got {!r}'.format(', '.join(SPIKEIN_QUANTIFIERS), SPIKEIN_QUANTIFIER))
if SPIKEIN_QUANTIFIER != 'skip' and not SPIKEIN_REF_SETS:
    raise ValueError('filter.spikein.quantifier is {} but filter.spikein.ref lists no spike-in sets'.format(SPIKEIN_QUANTIFIER))
if SPIKEIN_QUANTIFIER == 'skip' and SPIKEIN_REF_SETS:
    raise ValueError('filter.spikein.ref lists spike-in sets but filter.spikein.quantifier is skip; set bowtie (count and '
                     'remove them before quantification) or native (the quantifiers estimate them), or remove the sets')
if 'filter' in config['filter']['spikein']:
    raise ValueError('filter.spikein.filter is no longer used: filter.spikein.quantifier: bowtie counts and removes spike-in '
                     'reads before quantification, native leaves them to the quantifiers; remove filter.spikein.filter')


# A project's own FASTA files, by their name in the filter interim directory: rule spikein_project_fasta
# (filter/spikein_bowtie.smk) writes each in DNA letters there, so every route reads the same sequences.
SPIKEIN_CUSTOM = {'custom{}'.format(i): path for i, path in enumerate([n for n in SPIKEIN_REF_SETS if n not in SPIKEIN_SETS], 1)}


def native_spikein_fasta():
    """The spike-in FASTA files for quantifiers with their own spike-in support (UNITAS, miRge3): those of
    spikein_fasta() with filter.spikein.quantifier: native, none otherwise."""
    return spikein_fasta() if SPIKEIN_QUANTIFIER == 'native' else []


def spikein_fasta():
    """FASTA files of the configured spike-in sets, in order: gcfdb's for set names, the DNA copy for project files."""
    custom = {path: name for name, path in SPIKEIN_CUSTOM.items()}
    return [join(EXT_DIR, 'spikein', 'fasta', name + '.fa') if name in SPIKEIN_SETS
            else join(FILTER_INTERIM, 'spikein', 'sets', custom[name] + '.fa') for name in SPIKEIN_REF_SETS]
