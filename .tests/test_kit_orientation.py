"""rnaseq kit entries in libprep.config: read 1's orientation relative to the transcript (#362).

`read_orientation` sets salmon's library type and featureCounts' strand (rnaseq/rules/quant/salmon.smk,
featurecounts.smk): `reverse` gives `SR`/`ISR` and `-s 2`, `forward` gives `SF`/`ISF` and `-s 1`. `library_strand` says the
same in cDNA terms: `first strand` for `reverse`, `second strand` for `forward`.
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
STRAND = {'reverse': 'first strand', 'forward': 'second strand'}


def family(name):
    return re.sub(r' (SE|PE)$', '', name)


class RnaseqKitOrientationTests(unittest.TestCase):
    def test_every_rnaseq_kit_has_a_reviewed_orientation(self):
        self.assertEqual({family(name) for name in RNASEQ}, set(EXPECTED))

    def test_read_orientation_follows_the_kit_chemistry(self):
        for name, kit in RNASEQ.items():
            with self.subTest(kit=name):
                self.assertEqual(kit['read_orientation'], EXPECTED[family(name)])

    def test_library_strand_agrees_with_read_orientation(self):
        for name, kit in RNASEQ.items():
            with self.subTest(kit=name):
                self.assertEqual(kit['library_strand'], STRAND[kit['read_orientation']])


if __name__ == '__main__':
    unittest.main()
