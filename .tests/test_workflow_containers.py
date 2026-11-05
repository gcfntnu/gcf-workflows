"""No workflow-level container directive (#285).

A top-level `container:` applies to every later rule without its own container, across includes,
so the image a rule gets would depend on include order. Rules set their container explicitly.
"""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
SKIP = {'.dev', '.git', '.snakemake', '.tests'}
GLOBAL = re.compile(r'^(?:container|singularity)\s*:', re.M)


class WorkflowContainerTests(unittest.TestCase):
    def test_no_workflow_level_container(self):
        offenders = []
        for path in sorted(ROOT.rglob('*')):
            if (path.is_file() and (path.suffix in ('.smk', '.db') or path.name == 'Snakefile')
                    and not SKIP & set(path.relative_to(ROOT).parts) and GLOBAL.search(path.read_text(errors='replace'))):
                offenders.append(str(path.relative_to(ROOT)))
        self.assertEqual(offenders, [])

    def test_pattern_ignores_rule_directives(self):
        self.assertIsNone(GLOBAL.search("rule a:\n    container:\n        'docker://x:1'\n"))
        self.assertIsNotNone(GLOBAL.search("container: 'docker://x:1'\n"))


if __name__ == '__main__':
    unittest.main()
