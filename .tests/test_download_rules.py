"""Rules that download with wget or curl run without a container (#281).

They run in the BFQ image, which provides wget and curl. Default container images before
gcfntnu/python-default:3.12-1 have neither, so a `container:` on such a rule failed there with exit status 127
as soon as the download actually ran.
"""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
SKIP = {'.dev', '.git', '.snakemake', '.tests'}
RULE = re.compile(r'^([ \t]*)(?:rule|checkpoint)\s+(\w+)\s*:', re.M)
DOWNLOAD = re.compile(r'\b(wget|curl)\b')


def rule_files():
    for path in sorted(ROOT.rglob('*')):
        if path.is_file() and (path.suffix in ('.smk', '.db') or path.name == 'Snakefile') and not SKIP & set(path.relative_to(ROOT).parts):
            yield path


def rules(text):
    """(name, body) for every rule, including rules indented inside if blocks."""
    starts = list(RULE.finditer(text))
    for i, match in enumerate(starts):
        end = starts[i + 1].start() if i + 1 < len(starts) else len(text)
        indent = match.group(1)
        body = re.split(r'^' + re.escape(indent) + r'\S', text[match.end():end], flags=re.M)[0]
        yield match.group(2), indent, body


def directive(body, indent, name):
    return re.search(r'^' + re.escape(indent) + r'[ \t]+' + name + r'\s*:', body, re.M) is not None


class DownloadRuleTests(unittest.TestCase):
    def test_rules_that_call_wget_or_curl_have_no_container(self):
        downloads, offenders = 0, []
        for path in rule_files():
            for name, indent, body in rules(path.read_text(errors='replace')):
                command = body.split('shell:', 1)[-1] if 'shell:' in body else body.split('run:', 1)[-1] if 'run:' in body else ''
                if not DOWNLOAD.search(command):
                    continue
                downloads += 1
                if directive(body, indent, 'container') or directive(body, indent, 'singularity'):
                    offenders.append('{}: {}'.format(path.relative_to(ROOT), name))
        self.assertGreater(downloads, 50, 'the scan should find the reference getters')
        self.assertEqual(offenders, [], 'rules that download with wget/curl must not set a container')

    def test_scan_sees_container_and_indented_rules(self):
        text = ("rule plain:\n    output: 'a'\n    shell: 'wget -O {output} URL'\n"
                "if True:\n    rule nested:\n        output: 'b'\n        container:\n            'docker://example/x:1'\n"
                "        shell:\n            'curl -o {output} URL'\n")
        found = {name: directive(body, indent, 'container') for name, indent, body in rules(text)}
        self.assertEqual(found, {'plain': False, 'nested': True})


if __name__ == '__main__':
    unittest.main()
