"""UNITAS misc_RNA reads by gene family (#379), from the target names in unitas.hits_per_target.txt.

UNITAS reports Ensembl's misc_RNA biotype as one class; in blood-derived libraries it is mostly Y RNA fragments.
A misc_RNA target belongs to the first family whose rule matches the start of its name (the reference's Ensembl gene
name, e.g. `RNY4-ENSG00000252316.1`), ignoring case, so mouse and rat symbols (`Rny1`, `Rn7sk`, `Rn7sl1`) match too:

    Y_RNA  RNY..., Y_RNA         7SK  RN7SK..., 7SK         7SL  RN7SL..., Metazoa_SRP

Names that match no rule are `other`. A read UNITAS splits between targets keeps its fractional counts, so the
families add up to the misc_RNA class of the sample's annotation summary; check_total stops on a difference.
Imported by unitas_annotations.py and unitas_annotation2json.py.
"""
import csv
import re

FAMILIES = (('Y_RNA', re.compile(r'RNY|Y_RNA', re.I)), ('7SK', re.compile(r'RN7SK|7SK', re.I)),
            ('7SL', re.compile(r'RN7SL|Metazoa_SRP', re.I)))
OTHER = 'other'
NAMES = tuple(name for name, _ in FAMILIES) + (OTHER,)
HEADER = ['TRANSCRIPT_CLASS', 'TRANSCRIPT_NAME', 'READ_COUNT']


def family(target):
    """The misc_RNA family of a UNITAS target name."""
    return next((name for name, rule in FAMILIES if rule.match(target)), OTHER)


def misc_rna_families(path):
    """misc_RNA reads of one hits_per_target table by family, every family present."""
    counts = dict.fromkeys(NAMES, 0.0)
    with open(path, newline='') as fh:
        rows = csv.reader(fh, delimiter='\t')
        if next(rows, None) != HEADER:
            raise ValueError('{}: not a UNITAS hits_per_target table (expected columns {})'.format(path, ' '.join(HEADER)))
        for row in rows:
            if len(row) != 3:
                raise ValueError('{}: expected 3 columns, got {!r}'.format(path, row))
            if row[0] == 'misc_RNA':
                counts[family(row[1])] += float(row[2])
    return counts


def check_total(counts, misc_rna, path):
    """Fail unless the families add up to the misc_RNA class count of the annotation summary (0 when it has none)."""
    total = sum(counts.values())
    if abs(total - misc_rna) > 0.01 + 1e-6 * misc_rna:
        raise ValueError('{}: misc_RNA targets add up to {:.6g} reads, the annotation summary has {:.6g}'.format(path, total, misc_rna))

