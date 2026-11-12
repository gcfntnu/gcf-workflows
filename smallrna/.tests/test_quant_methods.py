"""Dry-run contract for the small-RNA `quant.method` selection of `quant_all`."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def configure(tmp_path, method, update=None):
    if not (tmp_path / 'Snakefile').exists():
        shutil.copytree(ROOT / 'smallrna/.tests', tmp_path, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
        (tmp_path / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                          f'include: {str(ROOT / "smallrna/smallrna.smk")!r}\n')
    settings = yaml.safe_load((tmp_path / 'config.yaml').read_text())
    settings.setdefault('quant', {})['method'] = method
    if update:
        update(settings)
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(settings))


def snakemake(tmp_path, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '--cores', '1', *args],
                            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
    return result.returncode, result.stdout + result.stderr


def dry_run(tmp_path, method, update=None, printshell=False, target='quant_all'):
    configure(tmp_path, method, update)
    return snakemake(tmp_path, '-n', *(['-p'] if printshell else []), target)


def record_as_done(tmp_path, *files):
    """Create files in order (inputs first) and record them as made by the current configuration."""
    for name in files:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.mkdir() if name.endswith('/') else path.touch()
        time.sleep(0.01)
    returncode, output = snakemake(tmp_path, '--touch', '--force', files[-1].rstrip('/'))
    assert returncode == 0, output


@pytest.mark.parametrize('method, rule', [('unitas', 'unitas_mirtable'), ('mirge3', 'mirge3_quant'),
                                          ('unitas,mirge3', 'mirge3_quant'), ('mirflowz', 'mirflowz_quant'),
                                          ('unitas,mirflowz', 'mirflowz_quant')])
def test_supported_methods_plan_their_quantifier(tmp_path, method, rule):
    returncode, output = dry_run(tmp_path, method)
    assert returncode == 0, output
    assert f'rule {rule}:' in output


@pytest.mark.parametrize('method', ['mirge', 'quickmirseq', 'seqbuster', 'unitass'])
def test_removed_or_unknown_methods_fail(tmp_path, method):
    returncode, output = dry_run(tmp_path, method)
    assert returncode != 0
    assert f'quant.method: unsupported value(s) {method}' in output


def mirge3_command(output):
    return next(line for line in output.splitlines() if 'miRge3.0 ' in line)


@pytest.mark.parametrize('kit, min_length', [('Bioo Scientific NEXTflex Small RNA-Seq Kit v3', 16),
                                             ('Bioo Scientific NEXTflex Small RNA-Seq Kit v4', 15),
                                             ('TruSeq Small RNA Library Prep Kit', 15)])
def test_mirge3_reads_the_filtered_reads_without_its_own_preprocessing(tmp_path, kit, min_length):
    """Like UNITAS and mirflowz: no adapter trimming, NextSeq trimming or UMI collapsing; the kit's minimum length."""
    update = lambda s: s.update(libprepkit=kit, machine='NextSeq 500')
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    quant = rule_block(output, 'mirge3_quant')
    assert all('data/tmp/smallrna/filter/fastp/{}_R1.fastq'.format(s) in quant.split('output:', 1)[0]
               for s in ('SR01', 'SR02', 'SR03', 'SR04'))
    assert 'ln -sf "$PWD"/data/tmp/smallrna/filter/fastp/SR01_R1.fastq data/tmp/smallrna/quant/mirge3/fastq/SR01.fastq' in quant
    command = ' '.join(mirge3_command(output).split())
    assert '--samples data/tmp/smallrna/quant/mirge3/fastq/SR01.fastq,' in command
    assert '--minimum-length {} '.format(min_length) in command
    for option in ('--adapter', '--nextseq-trim', '--uniq-mol-ids'):
        assert option not in command
    assert 'rule mirge3_symlink:' not in output


def test_mirge3_minimum_length_is_1_without_a_length_filter(tmp_path):
    """With -L the filter keeps every non-empty read, so miRge3 does too."""
    args = '-L -q 15 -u 40 -n 5 -e 0'
    update = lambda s: s.setdefault('filter', {}).setdefault('trim', {}).setdefault('fastp', {}).update(args=args)
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    assert '--minimum-length 1 ' in ' '.join(mirge3_command(output).split())


def test_mirge3_reads_every_lane_through_the_filter(tmp_path):
    def two_lanes(settings):
        sample = next(iter(settings['samples']))
        first = settings['samples'][sample]['R1']
        settings['samples'][sample]['R1'] = first + ',' + first.replace('_R1', '_R1_lane2')
        lane2 = tmp_path / settings['fastq_dir'] / first.replace('_R1', '_R1_lane2')
        shutil.copy(tmp_path / settings['fastq_dir'] / first, lane2)
    returncode, output = dry_run(tmp_path, 'mirge3', two_lanes, printshell=True)
    assert returncode == 0, output
    merge = [line for line in output.splitlines() if '_R1_lane2' in line]
    assert merge, 'the second lane was not merged'
    assert 'data/tmp/smallrna/filter/fastp/SR01_R1.fastq' in rule_block(output, 'mirge3_quant').split('output:', 1)[0]


@pytest.mark.parametrize('method, rule', [('mirge3', 'mirge3_quant'), ('unitas', 'unitas')])
def test_with_trimming_skipped_quantifiers_read_the_merged_reads(tmp_path, method, rule):
    """filter.trim.quantifier: skip gives the merged reads unchanged, uncompressed for the quantifiers."""
    update = lambda s: s.setdefault('filter', {}).setdefault('trim', {}).update(quantifier='skip')
    returncode, output = dry_run(tmp_path, method, update, printshell=True)
    assert returncode == 0, output
    # Jobs are planned in no fixed order, so look at every job of a rule.
    blocks = lambda name: [b.split('\nrule ', 1)[0] for b in output.split('rule {}:'.format(name))[1:]]
    assert any('data/tmp/smallrna/filter/fastq/SR01_R1.fastq' in b.split('output:', 1)[0] for b in blocks(rule))
    assert any('gzip -dc data/tmp/smallrna/filter/fastq/SR01_R1.fastq.gz > data/tmp/smallrna/filter/fastq/SR01_R1.fastq' in b
               for b in blocks('smallrna_untrimmed_fastq'))
    assert 'rule smallrna_fastp:' not in output


def test_mirge3_image_and_spikein_come_from_configuration(tmp_path):
    returncode, output = dry_run(tmp_path, 'mirge3', printshell=True)
    assert returncode == 0, output
    assert '--spikeIn' not in mirge3_command(output)
    assert 'quay.io/biocontainers/mirge3' not in (ROOT / 'smallrna/rules/quant/mirge3.smk').read_text()


def rule_block(output, rule):
    return output.split(f'rule {rule}:', 1)[1].split('\nrule ', 1)[0]


@pytest.mark.parametrize('ref, fasta', [(None, []), ('ERCC92', ['ERCC92']), ('smallrna_calibrators', ['smallrna_calibrators']),
                                        ('smallrna_calibrators, ERCC92', ['smallrna_calibrators', 'ERCC92'])])
def test_mirge3_spikein_index_follows_the_configured_sets(tmp_path, ref, fasta):
    """The spike-in index is built from the gcfdb sets of filter.spikein.ref before miRge3 runs, or removed."""
    update = lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref=ref)
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    project = rule_block(output, 'mirge3_project_library')
    command = next(line for line in project.splitlines() if 'mirge3_project_library.py' in line)
    assert ('--library data/ext/mirge3/09-11-2026/mouse --output data/tmp/smallrna/quant/mirge3/lib/mouse '
            '--record data/tmp/smallrna/quant/mirge3/spikein.tsv') in command
    assert command.split('--record data/tmp/smallrna/quant/mirge3/spikein.tsv', 1)[1].split() == \
        [f'data/ext/spikein/fasta/{name}.fa' for name in fasta]
    assert 'data/tmp/smallrna/quant/mirge3/lib/mouse' in rule_block(output, 'mirge3_quant').split('output:', 1)[0]
    command = ' '.join(mirge3_command(output).split())
    assert '--libraries-path data/tmp/smallrna/quant/mirge3/lib --organism-name mouse' in command
    assert ('--spikeIn' in command) == bool(fasta)


def test_unknown_spikein_set_fails_at_load(tmp_path):
    update = lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref='ERCC')
    returncode, output = dry_run(tmp_path, 'mirge3', update)
    assert returncode != 0
    assert "filter.spikein.ref: unsupported value 'ERCC'" in output


@pytest.mark.parametrize('db, name', [('mirbase', 'miRBase'), ('mirgenedb', 'MirGeneDB'), ('MirGeneDB', 'MirGeneDB')])
def test_mirge3_library_release_and_database_selection(tmp_path, db, name):
    update = lambda s: s['quant'].setdefault('mirge3', {}).update(db=db)
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    library = rule_block(output, 'mirge3_library')
    assert '--release-url mirge3/09-11-2026 --release 09-11-2026 --species mus_musculus --organism-name mouse' in library
    assert '--output data/ext/mirge3/09-11-2026/verify/miRge3_lib_mus_musculus.tsv' in library
    reference = rule_block(output, 'mirge3_reference')
    assert (f'--library data/ext/mirge3/09-11-2026/mouse --mirna-db {db.lower()} --release 09-11-2026 '
            '--verification data/ext/mirge3/09-11-2026/verify/miRge3_lib_mus_musculus.tsv') in reference
    command = ' '.join(mirge3_command(output).split())
    assert f'--mir-DB {name} --libraries-path data/tmp/smallrna/quant/mirge3/lib --organism-name mouse' in command
    assert 'data/tmp/smallrna/quant/mirge3/reference.tsv' in rule_block(output, 'quant_all')


def test_explicit_mirge3_release_selects_its_own_directory(tmp_path):
    update = lambda s: s.setdefault('db', {}).setdefault('mirge3', {}).update(release='01-12-2026')
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    assert '--release-url mirge3/01-12-2026' in rule_block(output, 'mirge3_library')
    assert '--library data/ext/mirge3/01-12-2026/mouse ' in rule_block(output, 'mirge3_project_library')


def test_sourceforge_library_stays_at_its_old_place(tmp_path):
    update = lambda s: s.setdefault('db', {}).setdefault('mirge3', {}).update(release='sourceforge')
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    assert 'rule mirge3_library:' not in output
    assert 'data/ext/mirge3/mus_musculus/mouse/index.Libs/mouse_snorna.3.ebwt' in rule_block(output, 'mirge3_prebuild')
    reference = rule_block(output, 'mirge3_reference')
    assert '--library data/ext/mirge3/mus_musculus/mouse --mirna-db mirbase --release sourceforge  --output' in reference
    assert '--library data/ext/mirge3/mus_musculus/mouse --output' in rule_block(output, 'mirge3_project_library')


LIBRARY = 'data/ext/mirge3/09-11-2026/verify/miRge3_lib_mus_musculus.tsv'


def test_switching_the_mirna_database_reruns_the_reference_check(tmp_path):
    """reference.tsv names the database it checked, so a switch in an existing workdir checks and records it again."""
    configure(tmp_path, 'mirge3', lambda s: s['quant'].setdefault('mirge3', {}).update(db='mirbase'))
    record_as_done(tmp_path, LIBRARY, 'data/tmp/smallrna/quant/mirge3/reference.tsv')
    returncode, output = dry_run(tmp_path, 'mirge3', lambda s: s['quant']['mirge3'].update(db='mirgenedb'),
                                 target='data/tmp/smallrna/quant/mirge3/reference.tsv')
    assert returncode == 0, output
    assert "Params have changed since last execution" in rule_block(output, 'mirge3_reference')
    assert "before: 'mirbase' now: 'mirgenedb'" in rule_block(output, 'mirge3_reference')


def test_changing_the_spikein_sets_rebuilds_the_project_library(tmp_path):
    """Each project's spike-in index is part of its own library, rebuilt when its sets change."""
    spikein = lambda ref: lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref=ref)
    configure(tmp_path, 'mirge3', spikein('smallrna_calibrators'))
    record_as_done(tmp_path, LIBRARY, 'data/ext/spikein/fasta/smallrna_calibrators.fa', 'data/ext/spikein/fasta/ERCC92.fa',
                   'data/tmp/smallrna/quant/mirge3/lib/mouse/', 'data/tmp/smallrna/quant/mirge3/spikein.tsv')
    returncode, output = dry_run(tmp_path, 'mirge3', target='data/tmp/smallrna/quant/mirge3/lib/mouse')
    assert returncode == 0, output
    assert 'Nothing to be done' in output
    returncode, output = dry_run(tmp_path, 'mirge3', spikein('ERCC92'), target='data/tmp/smallrna/quant/mirge3/lib/mouse')
    assert returncode == 0, output
    assert 'rule mirge3_project_library:' in output


@pytest.mark.parametrize('update, message', [
    (lambda s: s['quant'].setdefault('mirge3', {}).update(db='mirbase22'), "quant.mirge3.db: mirbase or mirgenedb, got 'mirbase22'"),
    (lambda s: s.setdefault('db', {}).setdefault('mirge3', {}).update(release=None), 'db.mirge3.release is not set')])
def test_invalid_mirge3_selection_fails_at_load(tmp_path, update, message):
    returncode, output = dry_run(tmp_path, 'mirge3', update)
    assert returncode != 0
    assert message in output


@pytest.mark.parametrize('mirna_db, feature_info', [
    ('mirbase', 'data/ext/mirbase/mus_musculus/release-23.1/anno/genes.tsv'),
    ('mirgenedb', 'data/ext/unitas/06-11-2026-mirgenedb/anno/mirgenedb_genes_mus_musculus.tsv')])
def test_unitas_mirna_feature_info_follows_the_dump(tmp_path, mirna_db, feature_info):
    """A MirGeneDB dump names miRNAs as MirGeneDB does, so the AnnData annotation comes from the dump's GFF."""
    def select(settings):
        settings.setdefault('db', {}).setdefault('unitas', {})['mirna_db'] = mirna_db
    returncode, output = dry_run(tmp_path, 'unitas', update=select, target='unitas_mirna_anndata')
    assert returncode == 0, output
    rule = output.split('rule unitas_mirna_anndata:', 1)[1]
    assert feature_info in rule.split('output:', 1)[0], output


@pytest.mark.parametrize('mirna_db, annotation', [
    ('mirbase', 'data/ext/mirbase/mus_musculus/release-23.1/anno/genes.gff3'),
    ('mirgenedb', 'data/ext/unitas/06-11-2026-mirgenedb/UNITAS_refdump_mus_musculus/gcf_mirgenedb.gff3')])
def test_unitas_reference_is_checked_before_counts_meet_annotations(tmp_path, mirna_db, annotation):
    """The dump is verified once in EXT_DIR, then checked against the annotation the AnnData joins it to."""
    def select(settings):
        settings.setdefault('db', {}).setdefault('unitas', {})['mirna_db'] = mirna_db
    returncode, output = dry_run(tmp_path, 'unitas', update=select, target='unitas_mirna_anndata', printshell=True)
    assert returncode == 0, output
    release = '06-11-2026-mirgenedb' if mirna_db == 'mirgenedb' else '05-11-2026'
    assert f'--output data/ext/unitas/{release}/verify/UNITAS_refdump_mus_musculus.tsv' in output
    check = output.split('rule unitas_reference:', 1)[1].split('rule ', 1)[0]
    assert f'--mirna-db {mirna_db} --annotation {annotation} --release {release}' in check
    assert 'data/tmp/smallrna/quant/unitas/reference.tsv' in output.split('rule unitas_mirna_anndata:', 1)[1].split('output:', 1)[0]
