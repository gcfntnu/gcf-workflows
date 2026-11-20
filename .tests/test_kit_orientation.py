"""rnaseq kit entries in libprep.config: read 1's orientation relative to the transcript (#362).

`read_orientation` sets salmon's library type and featureCounts' strand (rnaseq/rules/quant/salmon.smk,
featurecounts.smk): `reverse` gives `SR`/`ISR` and `-s 2`, `forward` gives `SF`/`ISF` and `-s 1`, `unstranded` gives `U`/`IU`
and `-s 0`. `library_strand` says the same in cDNA terms: `first strand` for `reverse`, `second strand` for `forward`,
`unstranded` for `unstranded`.
"""
import re
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
KITS = yaml.safe_load((ROOT / 'libprep.config').read_text())
RNASEQ = {name: kit for name, kit in KITS.items() if kit.get('workflow') == 'rnaseq'}

# Read 1 orientation per rnaseq kit (SE and PE entries alike), with its evidence. A new rnaseq kit needs a line here.
EXPECTED = {
    # dUTP marks the second cDNA strand, so read 1 is antisense. For Illumina Stranded mRNA Prep, salmon on earlier
    # facility analyses found nearly all fragments SR/ISR; the TruSeq Stranded kits use the same chemistry.
    'Illumina Stranded mRNA Prep': 'reverse',
    'Illumina Stranded Total mRNA Prep with Ribo-Zero Plus': 'reverse',
    'Illumina TruSeq Stranded mRNA Library Prep': 'reverse',
    'Illumina TruSeq Stranded Total RNA Library Prep (Globin)': 'reverse',
    # Lexogen SENSE: salmon run as SR on earlier facility analyses (mRNA V2, single end) mapped and assigned 86-91% of
    # the reads, which a forward library would not. SENSE Total uses the same library chemistry.
    'Lexogen SENSE mRNA-Seq Library Prep Kit V2': 'reverse',
    'Lexogen SENSE Total RNA-Seq Library Prep Kit (w/RiboCop rRNA Depletion Kit V1.2)': 'reverse',
    # Lexogen CORALL: read 1 sense, as all its entries are configured; no facility data.
    'Lexogen CORALL mRNA-Seq Library Prep Kit': 'forward',
    'CORALL Total RNA-Seq Library Prep Kit (w/RiboCop rRNA Depletion Kit V1.2)': 'forward',
}
STRAND = {'reverse': 'first strand', 'forward': 'second strand', 'unstranded': 'unstranded'}


def family(name):
    return re.sub(r' (SE|PE)$', '', name)


def mismatches(kits, expected):
    """Entries whose read_orientation is not the reviewed one for their family, or whose library_strand disagrees with
    their read_orientation."""
    found = []
    for name, kit in sorted(kits.items()):
        orientation, reviewed = kit.get('read_orientation'), expected.get(family(name))
        if orientation != reviewed:
            found.append('{}: read_orientation {} (reviewed: {})'.format(name, orientation, reviewed))
        if kit.get('library_strand') != STRAND.get(orientation):
            found.append('{}: library_strand {} with read_orientation {}'.format(name, kit.get('library_strand'), orientation))
    return found


class RnaseqKitOrientationTests(unittest.TestCase):
    def test_every_rnaseq_kit_has_a_reviewed_orientation(self):
        self.assertEqual({family(name) for name in RNASEQ}, set(EXPECTED))

    def test_rnaseq_kits_follow_their_reviewed_orientation_and_strand(self):
        self.assertEqual(mismatches(RNASEQ, EXPECTED), [])

    def test_an_unstranded_kit_records_library_strand_unstranded(self):
        reviewed = {'Example Kit': 'unstranded'}
        self.assertEqual(mismatches({'Example Kit SE': {'read_orientation': 'unstranded', 'library_strand': 'unstranded'}}, reviewed), [])
        self.assertEqual(mismatches({'Example Kit SE': {'read_orientation': 'unstranded', 'library_strand': 'first strand'}}, reviewed),
                         ['Example Kit SE: library_strand first strand with read_orientation unstranded'])


if __name__ == '__main__':
    unittest.main()
