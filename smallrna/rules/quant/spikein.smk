#-*- mode:snakemake -*-
"""Spike-in sets of `filter.spikein.ref`, shared by the UNITAS and miRge3 routes.

Spike-ins are a library prep choice. The sequences come from gcfdb (`EXT_DIR/spikein/fasta/<set>.fa`), whatever the
quantifier: UNITAS adds them to its reference with `-refseq`, miRge3 gets a spike-in index built from them.
"""

SPIKEIN_SETS = ('ERCC92', 'smallrna_calibrators')

def spikein_sets():
    """The comma-separated spike-in sets in filter.spikein.ref, in the given order; empty when none."""
    ref = config['filter']['spikein'].get('ref')
    names = [name.strip() for name in str(ref or '').split(',') if name.strip()]
    unknown = [name for name in names if name not in SPIKEIN_SETS]
    if unknown:
        raise ValueError('filter.spikein.ref: unsupported value {!r}; use a comma-separated list of {}'.format(ref, ', '.join(SPIKEIN_SETS)))
    return names

SPIKEIN_REF_SETS = spikein_sets()

def spikein_fasta():
    """gcfdb FASTA files of the configured spike-in sets."""
    return [join(EXT_DIR, 'spikein', 'fasta', name + '.fa') for name in SPIKEIN_REF_SETS]
