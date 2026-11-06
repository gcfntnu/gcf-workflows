"""Every literal `src_gcf(...)` path in the workflow files exists (#303).

`src_gcf(path)` joins `path` to the directory of the Snakemake file that calls it. A missing script only fails
when its rule runs, so a renamed script whose rule was not updated goes unnoticed in dry runs.
"""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
SKIP = {'.dev', '.git', '.snakemake', '.tests'}
SRC_GCF = re.compile(r'''src_gcf\(\s*(['"])([^'"]+)\1\s*\)''')

# Pre-existing missing paths, tracked in #306. Remove entries as they are fixed.
KNOWN_MISSING = {
    ('singlecell/rules/quant/preprocess.smk', 'scripts/preprocess_integrate_harmony.py'),
    ('singlecell/rules/quant/preprocess.smk', 'scripts/preprocess_integrate_scvi.py'),
}


def missing_paths(root=ROOT):
    missing = set()
    for path in sorted(root.rglob('*')):
        if (path.is_file() and (path.suffix in ('.smk', '.db') or path.name == 'Snakefile')
                and not SKIP & set(path.relative_to(root).parts)):
            for match in SRC_GCF.finditer(path.read_text(errors='replace')):
                if not (path.parent / match.group(2)).exists():
                    missing.add((str(path.relative_to(root)), match.group(2)))
    return missing


class WorkflowScriptTests(unittest.TestCase):
    def test_src_gcf_paths_exist(self):
        self.assertEqual(missing_paths(), KNOWN_MISSING)

    def test_pattern(self):
        self.assertEqual(SRC_GCF.findall('script = src_gcf("scripts/a.py"),\nx = src_gcf(join("a", "b"))'), [('"', 'scripts/a.py')])


if __name__ == '__main__':
    unittest.main()
