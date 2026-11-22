"""Dry-run contract for the small-RNA `quant.method` selection of `quant_all`."""
import hashlib
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def mirbase_gff_id(root=ROOT):
    """MIRBASE_GFF_ID (gcfdb/mirbase.db, #438) of the checkout at `root`: a digest of the correction table and script, which
    names the corrected GFF3 and the shared files built from it (#440)."""
    digest = hashlib.sha256()
    for name in ('gcfdb/tables/mirbase_gff_corrections.tsv', 'gcfdb/scripts/mirbase_gff_correct.py'):
        digest.update(hashlib.sha256((Path(root) / name).read_bytes()).digest())
    return digest.hexdigest()[:12]


MIRBASE_GFF_ID = mirbase_gff_id()


def corrected_mirbase_gff(ext='data/ext', organism='mus_musculus', release='23.1'):
    """The corrected miRBase GFF3, named after the correction table and script (MIRBASE_GFF_ID in gcfdb/mirbase.db, #438)."""
    return f'{ext}/mirbase/{organism}/release-{release}/anno/genes.{MIRBASE_GFF_ID}.gff3'


MIRBASE_GFF = corrected_mirbase_gff()


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
                                          ('unitas,mirflowz', 'mirflowz_quant'), ('biocore_htseq', 'biocore_htseq_count'),
                                          ('unitas,biocore_htseq', 'biocore_htseq_count'),
                                          ('biocore_miraligner', 'biocore_miraligner_count'), ('unitas,biocore_miraligner', 'biocore_miraligner_count')])
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
                                             ('Bioo Scientific NEXTflex Small RNA-Seq Kit v4', 16),
                                             ('TruSeq Small RNA Library Prep Kit', 16)])
def test_mirge3_reads_the_filtered_reads_without_its_own_preprocessing(tmp_path, kit, min_length):
    """Like UNITAS and mirflowz: no adapter trimming, NextSeq trimming or UMI collapsing; the kit's minimum length."""
    update = lambda s: s.update(libprepkit=kit, machine='NextSeq 500')
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    quant = rule_block(output, 'mirge3_quant')
    assert all('data/tmp/smallrna/filter/cutadapt/{}_R1.fastq'.format(s) in quant.split('output:', 1)[0]
               for s in ('SR01', 'SR02', 'SR03', 'SR04'))
    assert 'ln -sf "$PWD"/data/tmp/smallrna/filter/cutadapt/SR01_R1.fastq data/tmp/smallrna/quant/mirge3/fastq/SR01.fastq' in quant
    command = ' '.join(mirge3_command(output).split())
    assert '--samples data/tmp/smallrna/quant/mirge3/fastq/SR01.fastq,' in command
    assert '--minimum-length {} '.format(min_length) in command
    for option in ('--adapter', '--nextseq-trim', '--uniq-mol-ids'):
        assert option not in command
    assert '--crThreshold 0 ' in command  # no hidden canonical-ratio filter (#347)
    assert 'rule mirge3_symlink:' not in output


def test_mirge3_minimum_length_is_1_without_a_length_filter(tmp_path):
    """With fastp's -L the filter keeps every non-empty read, so miRge3 does too."""
    update = lambda s: s.setdefault('filter', {}).update(trim={'quantifier': 'fastp', 'fastp': {'args': '-L -q 15 -u 40 -n 5 -e 0'}})
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
    assert 'data/tmp/smallrna/filter/cutadapt/SR01_R1.fastq' in rule_block(output, 'mirge3_quant').split('output:', 1)[0]


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
    """With filter.spikein.quantifier: native, the spike-in index is built from the gcfdb sets of filter.spikein.ref
    before miRge3 runs; without sets there is none."""
    update = lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref=ref, **({'quantifier': 'native'} if ref else {}))
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    project = rule_block(output, 'mirge3_project_library')
    command = next(line for line in project.splitlines() if 'mirge3_project_library.py' in line)
    assert ('--library data/ext/mirge3/20-11-2026-2/mouse --output data/tmp/smallrna/quant/mirge3/lib/mouse '
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
    assert '--release-url mirge3/20-11-2026-2 --release 20-11-2026-2 --species mus_musculus --organism-name mouse' in library
    assert '--output data/ext/mirge3/20-11-2026-2/verify/miRge3_lib_mus_musculus.tsv' in library
    reference = rule_block(output, 'mirge3_reference')
    assert (f'--library data/ext/mirge3/20-11-2026-2/mouse --mirna-db {db.lower()} --release 20-11-2026-2 '
            '--verification data/ext/mirge3/20-11-2026-2/verify/miRge3_lib_mus_musculus.tsv') in reference
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


LIBRARY = 'data/ext/mirge3/20-11-2026-2/verify/miRge3_lib_mus_musculus.tsv'


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
    spikein = lambda ref: lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref=ref, quantifier='native')
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
    ('mirbase', f'data/ext/mirbase/mus_musculus/release-23.1/anno/genes.{MIRBASE_GFF_ID}.tsv'),
    ('mirgenedb', 'data/ext/unitas/20-11-2026-mirgenedb/anno/mirgenedb_genes_mus_musculus.tsv')])
def test_unitas_mirna_feature_info_follows_the_dump(tmp_path, mirna_db, feature_info):
    """A MirGeneDB dump names miRNAs as MirGeneDB does, so the AnnData annotation comes from the dump's GFF."""
    def select(settings):
        settings.setdefault('db', {}).setdefault('unitas', {})['mirna_db'] = mirna_db
    returncode, output = dry_run(tmp_path, 'unitas', update=select, target='unitas_mirna_anndata')
    assert returncode == 0, output
    rule = output.split('rule unitas_mirna_anndata:', 1)[1]
    assert feature_info in rule.split('output:', 1)[0], output


@pytest.mark.parametrize('mirna_db, annotation', [
    ('mirbase', MIRBASE_GFF),
    ('mirgenedb', 'data/ext/unitas/20-11-2026-mirgenedb/UNITAS_refdump_mus_musculus/gcf_mirgenedb.gff3')])
def test_unitas_reference_is_checked_before_counts_meet_annotations(tmp_path, mirna_db, annotation):
    """The dump is verified once in EXT_DIR, then checked against the annotation the AnnData joins it to."""
    def select(settings):
        settings.setdefault('db', {}).setdefault('unitas', {})['mirna_db'] = mirna_db
    returncode, output = dry_run(tmp_path, 'unitas', update=select, target='unitas_mirna_anndata', printshell=True)
    assert returncode == 0, output
    release = '20-11-2026-mirgenedb' if mirna_db == 'mirgenedb' else '20-11-2026'
    assert f'--output data/ext/unitas/{release}/verify/UNITAS_refdump_mus_musculus.tsv' in output
    check = output.split('rule unitas_reference:', 1)[1].split('rule ', 1)[0]
    assert f'--mirna-db {mirna_db} --annotation {annotation} --release {release}' in check
    assert 'data/tmp/smallrna/quant/unitas/reference.tsv' in output.split('rule unitas_mirna_anndata:', 1)[1].split('output:', 1)[0]


@pytest.mark.parametrize('db, name', [('mirbase', 'miRBase'), ('mirgenedb', 'MirGeneDB')])
def test_mirge3_counts_are_made_from_its_mapped_reads(tmp_path, db, name):
    """miR.Counts.csv comes from mapped.csv with the library's merge file; miRge3's own table is kept (#349)."""
    returncode, output = dry_run(tmp_path, 'mirge3', lambda s: s['quant'].setdefault('mirge3', {}).update(db=db),
                                 printshell=True)
    assert returncode == 0, output
    quant = rule_block(output, 'mirge3_quant')
    assert 'data/tmp/smallrna/quant/mirge3/miR.Counts.mirge3.csv' in quant and 'quant/mirge3/mapped.csv' in quant
    counts = rule_block(output, 'mirge3_counts')
    assert 'output: data/tmp/smallrna/quant/mirge3/miR.Counts.csv' in counts
    assert f'--merges data/tmp/smallrna/quant/mirge3/lib/mouse/annotation.Libs/mouse_merges_{name}.csv' in counts


def test_recorded_routes_write_their_mirbase_release_and_table_checksum(tmp_path):
    """quant/<route>/reference.tsv for mirflowz, TEsmall, tiny-count, nf-core/smrnaseq, biocore_htseq and biocore_miraligner (#407),
    written after the route's count table with its SHA-256: db.mirbase.release, the mirflowz reference, gcfdb's miRBase
    files or miraligner's database written from them; NA for a db.mirflowz.custom reference; UNITAS and miRge3 write their
    own."""
    returncode, output = dry_run(tmp_path, 'unitas,mirflowz,tesmall,tinyrna,smrnaseq,biocore_htseq,biocore_miraligner', printshell=True)
    assert returncode == 0, output
    blocks = {re.search(r'output: data/tmp/smallrna/quant/(\w+)/reference.tsv', b).group(1): b for b in output.split('rule quant_reference:')[1:]}
    assert sorted(blocks) == ['biocore_htseq', 'biocore_miraligner', 'mirflowz', 'smrnaseq', 'tesmall', 'tinyrna']
    tables = {'mirflowz': 'mirflowz/all_mirna_counts.tab', 'tesmall': 'tesmall/tesmall_counts.tsv', 'tinyrna': 'tinyrna/tinyrna_counts.tsv',
              'smrnaseq': 'smrnaseq/mirna.tsv', 'biocore_htseq': 'biocore_htseq/htseq_counts.tsv', 'biocore_miraligner': 'biocore_miraligner/mirna.tsv'}
    for route, block in blocks.items():
        table = f'data/tmp/smallrna/quant/{tables[route]}'
        assert f'input: {table}' in block  # written after its count table
        assert f'sha=$(sha256sum {table} | cut' in block and f'"$sha" {table} > data/tmp/smallrna/quant/{route}/reference.tsv' in block
        reference = {'smrnaseq': 'data/ext/mirbase/mus_musculus/release-23.1',
                     'biocore_miraligner': f'data/ext/biocore_miraligner/mirbase-23.1/mus_musculus/gff-{MIRBASE_GFF_ID}/db.v1'}.get(
            route, 'data/ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus')
        # and the identity of the GFF its annotation was built from (#440), which names its gff-<ID> files
        assert f'mirbase db.mirbase.release 23.1 db.mirbase.release {reference} {MIRBASE_GFF_ID} "$sha"' in block
        assert "reference\\t%s\\tconfiguration\\nmirna_gff_id\\t%s\\tconfiguration" in block
    refs = tmp_path / 'refs'
    refs.mkdir()
    files = {key: refs / name for key, name in (('genome', 'genome.fa.gz'), ('gtf', 'genes.gtf.gz'), ('mirna_gff', 'mirna.gff3'))}
    for path in files.values():
        path.touch()
    custom = lambda s: s.setdefault('db', {}).setdefault('mirflowz', {}).update(custom=dict(name='reduced', **{k: str(v) for k, v in files.items()}))
    returncode, output = dry_run(tmp_path, 'tinyrna', custom, printshell=True)
    assert returncode == 0, output
    assert 'NA db.mirflowz.custom NA db.mirflowz.custom data/ext/mirflowz/custom-reduced/mus_musculus custom "$sha"' in rule_block(output, 'quant_reference')
    returncode, output = dry_run(tmp_path, 'unitas,mirge3')
    assert returncode == 0, output
    assert 'quant_reference' not in output


def test_route_steps_write_the_benchmark_files_the_benchmark_reads(tmp_path):
    """Every step benchmark/resources.py counts for a method declares its Snakemake benchmark file where resources.jobs
    reads it (#386): quant/benchmarks/<rule>/<sample>.tsv per sample, <rule>.tsv per project; the shared reference builds
    in their reference's benchmarks/ folder."""
    sys.path.insert(0, str(ROOT / 'smallrna/benchmark'))
    import resources
    returncode, output = dry_run(tmp_path, ','.join(resources.STEPS))
    assert returncode == 0, output
    declared = {}
    for block in output.split('\nrule ')[1:]:
        match = re.search(r'^    benchmark: (\S+)$', block, re.M)
        if match:
            declared.setdefault(block.split(':', 1)[0], set()).add(match.group(1))
    base = 'data/tmp/smallrna/quant/benchmarks/'
    samples = {Path(p).stem for p in declared['unitas']}
    assert len(samples) == 4
    for method, steps in resources.STEPS.items():
        for step, per_sample, _ in steps:
            expected = {f'{base}{step}/{s}.tsv' for s in samples} if per_sample else {f'{base}{step}.tsv'}
            assert declared.get(step) == expected, (method, step)
    # shared builds: in the reference they built, so a file never describes another reference
    reference = 'data/ext/{}/mirbase-23.1_ensembl-116/mus_musculus'
    assert declared['mirflowz_prepare'] == {reference.format('mirflowz') + '/benchmarks/mirflowz_prepare.tsv'}
    assert declared['mirflowz_annotation'] == {reference.format('mirflowz') + f'/gff-{MIRBASE_GFF_ID}/benchmarks/mirflowz_annotation.tsv'}
    assert declared['tesmall_genome'] == {reference.format('tesmall') + '/benchmarks/tesmall_genome.tsv'}
    assert declared['tesmall_db'] == {reference.format('tesmall') + f'/gff-{MIRBASE_GFF_ID}/unitas-20-11-2026/benchmarks/tesmall_db.tsv'}
    assert declared['biocore_htseq_bowtie2_index'] == {reference.format('biocore_htseq') + '/benchmarks/biocore_htseq_bowtie2_index.tsv'}
    assert declared['biocore_miraligner_db'] == {
        f'data/ext/biocore_miraligner/mirbase-23.1/mus_musculus/gff-{MIRBASE_GFF_ID}/db.v1/benchmarks/biocore_miraligner_db.tsv'}


ALL_ROUTES = 'unitas,mirge3,mirflowz,tesmall,tinyrna,smrnaseq,biocore_htseq,biocore_miraligner'


def settings_of(output):
    """{route: (settings [(field, value, source)], the quant_settings command's words)} of a dry run printed with -p."""
    import json
    import shlex
    out = {}
    for block in output.split('rule quant_settings:')[1:]:
        route = re.search(r'output: data/tmp/smallrna/quant/(\w+)/settings.tsv', block).group(1)
        line = next(line for line in block.splitlines() if 'quant_settings.py' in line)
        words = shlex.split(line.split('Shell command:', 1)[-1])
        out[route] = ([tuple(row) for row in json.loads(words[words.index('--settings') + 1])], words)
    return out


def test_each_route_records_the_settings_its_commands_use(tmp_path):
    """#417: quant_settings writes quant/<route>/settings.tsv for every selected route, after its count table, from the
    same values its commands use, so the record and the run cannot differ."""
    returncode, output = dry_run(tmp_path, ALL_ROUTES, printshell=True)
    assert returncode == 0, output
    settings = settings_of(output)
    assert sorted(settings) == sorted(ALL_ROUTES.split(','))
    fields = {route: {f: (v, s) for f, v, s in rows} for route, (rows, _) in settings.items()}
    assert fields['tinyrna']['mismatch'] == ('0', 'quant.tinyrna.mismatch') and fields['tinyrna']['normalize_by_genomic_hits'][0] == 'True'
    assert fields['tesmall']['order'] == ('miRNA,hairpin,structural_RNA,exon,TE,intron,piRNA_cluster', 'quant.tesmall.order')
    assert '--crThreshold 0' in fields['mirge3']['arguments'][0] and fields['mirge3']['mir_db'] == ('miRBase', 'quant.mirge3.db')
    assert fields['unitas']['refseq'][0] == 'none'
    # the record holds what the commands run with
    for route, rule, field in (('tinyrna', 'tinyrna_align', 'bowtie'), ('smrnaseq', 'smrnaseq_align', 'bowtie'),
                               ('tesmall', 'tesmall_quant', 'arguments'), ('biocore_htseq', 'biocore_htseq_annotate', 'htseq_count'),
                               ('biocore_miraligner', 'biocore_miraligner_align', 'miraligner')):
        recorded = fields[route][field][0].split(' (')[0].replace('--no-trim ', '').replace(' -freq', '')
        assert recorded in ' '.join(rule_block(output, rule).split()), (route, recorded)
    mirflowz_words = settings['mirflowz'][1]
    assert mirflowz_words[mirflowz_words.index('--recorded') + 1].endswith('mirflowz_config.yaml:minimum_length,overlap,max_n,max_length_reads,nh,'
                                                                          'extension,error_rate,q_value,p_value')
    sys.path.insert(0, str(ROOT / 'smallrna/rules/quant/scripts'))
    import mirge3_counts
    assert 'needs {} exact reads'.format(mirge3_counts.MIN_PROJECT_EXACT) in fields['mirge3']['project_rule'][0]


def test_quant_settings_script_binds_the_record_to_its_table(tmp_path):
    import json
    script = ROOT / 'smallrna/rules/quant/scripts/quant_settings.py'
    table = tmp_path / 'all_mirna_counts.tab'
    table.write_text('ID\tS1\nmiR-1\t5\n')
    run_config = tmp_path / 'mirflowz_config.yaml'
    run_config.write_text('nh: 100\nmax_length_reads: 30\nerror_rate: 0.1\nmir_list: ["isomir", "mirna"]\n')
    out = tmp_path / 'settings.tsv'
    cmd = [sys.executable, str(script), '--table', str(table), '--settings', json.dumps([['adapter', 'TGGAATTC', 'kit']]),
           '--recorded', '{}:nh,max_length_reads,error_rate'.format(run_config), '--output', str(out)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert out.read_text() == ('field\tvalue\tsource\nadapter\tTGGAATTC\tkit\nnh\t100\tmirflowz_config.yaml\n'
                               'max_length_reads\t30\tmirflowz_config.yaml\nerror_rate\t0.1\tmirflowz_config.yaml\n'
                               'table_sha256\t{}\t{}\n'.format(hashlib.sha256(table.read_bytes()).hexdigest(), table))
    for settings, recorded, message in ((json.dumps([['a', '1', 'x'], ['a', '2', 'x']]), None, 'fields must be unique'),
                                        (json.dumps([['a', 'x\ty', 's']]), None, 'a tab or a line break'),
                                        (json.dumps([['a', '1', 'x']]), '{}:nh,absent'.format(run_config), 'has no absent')):
        result = subprocess.run(cmd[:4] + ['--settings', settings] + (['--recorded', recorded] if recorded else []) + ['--output', str(out)],
                                capture_output=True, text=True)
        assert result.returncode != 0 and message in result.stderr, result.stderr
