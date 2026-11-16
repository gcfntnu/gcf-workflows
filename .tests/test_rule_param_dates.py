"""Rule params must not depend on the current date (#330).

Snakemake reruns a finished job when its params change, so a param such as datetime.now().strftime(...) made completed
references and everything downstream rerun when a workflow was resumed on a later day. Reference rules record the
download date when they run, in the shell command (`$(date +%d-%m-%Y)`), instead.
"""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
CLOCK = re.compile(r'datetime\.now\(|datetime\.today\(|date\.today\(|time\.time\(|strftime\(')
RULE_FILES = [p for p in ROOT.rglob('*') if p.suffix in ('.smk', '.db') or p.name == 'Snakefile']
RULE_FILES = [p for p in RULE_FILES if not any(part.startswith('.') and part not in ('.tests',) for part in p.relative_to(ROOT).parts)]


def params_lines(text):
    """(line number, line) of every line inside a rule's params: block."""
    inside, indent = False, 0
    for number, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        current = len(line) - len(line.lstrip())
        if inside and stripped and current <= indent:
            inside = False
        if re.match(r'^\s+params:\s*(#.*)?$', line):
            inside, indent = True, current
            continue
        if inside:
            yield number, line


class RuleParamDates(unittest.TestCase):
    def test_no_rule_param_reads_the_clock(self):
        self.assertTrue(any(p.name == 'mirbase.db' for p in RULE_FILES))
        found = ['{}:{}: {}'.format(p.relative_to(ROOT), n, line.strip())
                 for p in RULE_FILES for n, line in params_lines(p.read_text(errors='replace')) if CLOCK.search(line)]
        self.assertEqual(found, [])

    def test_params_lines_finds_a_date_param(self):
        text = 'rule a:\n    params:\n        url = "u",\n        date = datetime.now().strftime("%d")\n    shell:\n        "x"\n'
        self.assertEqual([n for n, line in params_lines(text) if CLOCK.search(line)], [4])


SNAKEFILE = textwrap.dedent('''
    rule download:
        output: "ref.txt"
        params:
            url = config["url"],
            **({"date": config["date"]} if "date" in config else {})
        log: "ref.log"
        shell: 'echo {params.url} > {output} && echo "ref,{params.url},$(date +%d-%m-%Y)" > {log}'
''')


class RerunOnLaterDay(unittest.TestCase):
    """With Snakemake's default rerun triggers, a date param reruns a finished job; a date written at run time does not,
    and a meaningful param change still does."""

    def snakemake(self, work, *config, dry=True):
        command = [str(Path(sys.executable).with_name('snakemake')), '--cores', '1', '--config', *config]
        result = subprocess.run(command + (['-n'] if dry else []), cwd=work, capture_output=True, text=True, timeout=120,
                                env={k: v for k, v in os.environ.items() if not k.startswith('GCF_')})
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout + result.stderr

    def test_rerun_triggers(self):
        with tempfile.TemporaryDirectory() as work:
            Path(work, 'Snakefile').write_text(SNAKEFILE)
            self.snakemake(work, 'url=https://example.org/a', dry=False)
            self.assertRegex(Path(work, 'ref.log').read_text(), r'^ref,https://example.org/a,\d\d-\d\d-\d{4}\n$')  # date of the run
            self.assertIn('Nothing to be done', self.snakemake(work, 'url=https://example.org/a'))
            self.assertIn('rule download', self.snakemake(work, 'url=https://example.org/b'))

        with tempfile.TemporaryDirectory() as work:  # the former pattern: the date as a param
            Path(work, 'Snakefile').write_text(SNAKEFILE)
            self.snakemake(work, 'url=https://example.org/a', 'date=16-11-2026', dry=False)
            self.assertIn('Nothing to be done', self.snakemake(work, 'url=https://example.org/a', 'date=16-11-2026'))
            self.assertIn('rule download', self.snakemake(work, 'url=https://example.org/a', 'date=17-11-2026'))


if __name__ == '__main__':
    unittest.main()
