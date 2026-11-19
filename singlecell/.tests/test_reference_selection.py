"""Plan the BFQ reporting target and inspect the tools' reference arguments with an empty cache."""
import hashlib
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

import pytest
import yaml

from test_velocity import TENX_KIT, WORKFLOW, protected_env

DB_DEFAULTS = yaml.safe_load((WORKFLOW.parents[1] / 'main.config').read_text())['db']
COUNT_TARGET = '_tmp_libA/outs/web_summary.html'
MOCK_TOOLS = WORKFLOW.parents[1] / 'gcfdb/.tests/mock_tools'


def dry_run(work, organism, target='multiqc_report', expect_success=True, workflow_file=WORKFLOW, kit=TENX_KIT, **settings):
    config = dict(workflow=workflow_file.stem, organism=organism, project_id=['SYNTHETIC'], libprepkit=kit,
                  machine='NextSeq 500', read_geometry=[28, 90],
                  samples={'libA': {'R1': 'libA_R1.fastq.gz', 'R2': 'libA_R2.fastq.gz'}}, **settings)
    (work / 'config.yaml').write_text(yaml.safe_dump(config))
    (work / 'pep.yaml').write_text('pep_version: 2.0.0\nsample_table: samples.csv\n')
    (work / 'samples.csv').write_text('sample_name,R1,R2,paired_end\nlibA,libA_R1.fastq.gz,libA_R2.fastq.gz,1\n')
    fastqs = work / 'data/raw/fastq'
    fastqs.mkdir(parents=True)
    for read in config['samples']['libA'].values():
        (fastqs / read).touch()
    (work / 'Snakefile').write_text(f'pepfile: "pep.yaml"\nconfigfile: "config.yaml"\ninclude: {str(workflow_file)!r}\n')
    assert not (work / 'data/ext').exists()
    return plan_existing(work, target, expect_success)


def plan_existing(work, target, expect_success=True):
    result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', target],
                            cwd=work, env=protected_env(work), capture_output=True, text=True, timeout=120)
    log = result.stdout + result.stderr
    assert (result.returncode == 0) == expect_success, log
    return log


def arguments(log, command, flag):
    """Read only the argument values actually passed in the expanded shell commands."""
    values = set()
    for line in log.splitlines():
        if line.startswith('Shell command:') and command in line:
            args = shlex.split(line)
            values.add(args[args.index(flag) + 1])
    assert values, f'Missing {command} {flag} in dry run:\n{log}'
    return values


def assert_build_and_count(log, ref):
    assert arguments(log, 'cellranger mkref ', '--fasta') == {f'{ref}/fasta/genome.fa'}
    assert arguments(log, 'cellranger mkref ', '--genes') == {f'{ref}/anno/genes.gtf'}
    assert arguments(log, 'cellranger mkref ', '--cache') == {f'{ref}/cellranger/builds'}
    receipt, = arguments(log, 'cellranger mkref ', '--receipt')
    expected = Path(receipt).parent / 'index'  # A managed build pins the completed generation here.
    if '--expected-source' in next(line for line in log.splitlines() if 'cellranger mkref ' in line):
        source, = arguments(log, 'cellranger mkref ', '--expected-source')
        build_key = Path(receipt).parent.name
        expected = Path(ref) / 'cellranger/builds' / source / build_key
    assert arguments(log, 'cellranger count ', '--transcriptome') == {str(expected)}


def assert_gene_metadata(log, ref):
    assert arguments(log, 'bedtools getfasta ', '-fi') == {f'{ref}/fasta/genome.fa'}
    assert arguments(log, 'bedtools getfasta ', '-bed') == {f'{ref}/anno/exons.bed'}
    assert arguments(log, 'gtf2tsv.py ', '--fasta') == {f'{ref}/fasta/exons.fa'}
    assert arguments(log, 'gtf2tsv.py ', '--gtf') == {f'{ref}/anno/genes.gtf'}
    assert arguments(log, 'gtf2tsv.py ', '--output') == {f'{ref}/anno/genes.tsv'}
    assert arguments(log, 'convert_scanpy.py ', '--feature-info') == {f'{ref}/anno/genes.tsv'}


@pytest.mark.parametrize('organism,assembly,settings', [
    ('salmo_salar', 'Ssal_v3.1', {}),
    ('danio_rerio', 'GRCz11', {'db': {'10xgenomics': {'release': None}}}),
])
def test_ensembl_fallback_through_reporting(tmp_path, organism, assembly, settings):
    log = dry_run(tmp_path, organism, **settings)
    release = DB_DEFAULTS['ensembl']['release']
    ref = f'data/ext/ensembl/release-{release}/{organism}/{assembly}'
    assert 'rule multiqc_report:' in log
    assert f'source=ensembl release={release} assembly={assembly}' in log
    assert f'{organism.capitalize()}.{assembly}.dna.toplevel.fa.gz' in log
    assert f'{organism.capitalize()}.{assembly}.{release}.gtf.gz' in log
    assert_build_and_count(log, ref)
    # Keep the assembly in the paths and provenance, but pass a valid Cell Ranger genome identifier.
    assert arguments(log, 'cellranger mkref ', '--genome') == {assembly.replace('.', '_')}
    receipt, = arguments(log, 'cellranger mkref ', '--receipt')
    assert_managed_gene_metadata(log, Path(receipt))
    assert 'rule txgenomics_org_prebuild:' not in log
    assert 'refdata-gex-' not in log


def test_ensembl_genomes_fallback_and_workdir_release_change(tmp_path):
    organism, assembly = 'arabidopsis_thaliana', 'TAIR10'
    release = DB_DEFAULTS['ensemblgenomes']['release']
    log = dry_run(tmp_path, organism)
    for selected in (release, str(int(release) - 1)):
        if selected != release:
            configfile = tmp_path / 'config.yaml'
            config = yaml.safe_load(configfile.read_text())
            config['db'] = {'ensemblgenomes': {'release': selected}}
            configfile.write_text(yaml.safe_dump(config))
            log = plan_existing(tmp_path, 'multiqc_report')

        archive = f'https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-{selected}/plants'
        genome_url = f'{archive}/fasta/{organism}/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'
        gtf_url = f'{archive}/gtf/{organism}/Arabidopsis_thaliana.TAIR10.{selected}.gtf.gz'
        assert f'-O- {genome_url} |' in log and f'-O- {gtf_url} |' in log
        assert f'Ensembl DNA,release-{selected},{genome_url}' in log and f'Ensembl GTF,release-{selected},{gtf_url}' in log
        assert f'source=ensemblgenomes release={selected} assembly={assembly}' in log and 'division=plants' in log
        location = f'ensemblgenomes/release-{selected}/plants/{organism}/{assembly}'
        assert_build_and_count(log, f'data/ext/{location}')
        receipt, = arguments(log, 'cellranger mkref ', '--receipt')
        assert f'/{location}/' in receipt
        assert_managed_gene_metadata(log, Path(receipt))
        assert 'rule multiqc_report:' in log and 'refdata-gex-' not in log


def shell_command(log, program):
    lines = [line.removeprefix('Shell command: ') for line in log.splitlines() if line.startswith('Shell command: ') and program in line]
    assert len(lines) == 1, lines
    return shlex.split(lines[0])


def assert_managed_gene_metadata(log, receipt):
    view = receipt.parent
    assert_gene_metadata(log, view / 'annotation')
    stage = shell_command(log, 'stage_project_annotation.sh')
    assert stage[2:] == [str(view / 'index/gcf-source/fasta/genome.fa'), str(view / 'index/gcf-source/anno/genes.gtf'),
                         str(view / 'annotation'), '--manifest', str(receipt)]


@pytest.mark.parametrize('organism,provider_config,reference,changed,removed', [
    ('salmo_salar', {}, f"ensembl/release-{DB_DEFAULTS['ensembl']['release']}/salmo_salar/Ssal_v3.1", 'anno/genes.gtf', '.'),
    ('rangifer_tarandus', {'db': {'reference_db': 'custom_reference', 'custom_reference': {'assembly': 'Rtarandus_2.1'}}},
     'custom_reference/rangifer_tarandus/Rtarandus_2.1', 'fasta/genome.fa', 'gcf-mkref.json'),
])
def test_managed_count_and_gene_metadata_share_verified_snapshot(tmp_path, organism, provider_config, reference, changed, removed):
    shared = tmp_path / 'ext'
    work = tmp_path / 'work'
    work.mkdir()
    cold = dry_run(work, organism, ext_dir=str(shared), **provider_config)
    assert 'rule multiqc_report:' in cold
    receipt = Path(next(iter(arguments(cold, 'cellranger mkref ', '--receipt'))))
    assert_managed_gene_metadata(cold, receipt)

    # Stand in for the provider download and mkref binary, then run the production publisher and copy helper.
    source = shared / reference
    write_panel(tmp_path, source)
    for path in ('fasta/genome.fa', 'anno/genes.gtf'):
        os.utime(source / path, (100, 100))
    env = protected_env(work)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    publish = shell_command(cold, 'cellranger mkref ')
    subprocess.run([sys.executable, *publish[1:]], cwd=work, env=env, check=True)

    warm = plan_existing(work, 'multiqc_report')
    assert 'rule cellranger_mkref:' not in warm
    assert arguments(warm, 'cellranger count ', '--transcriptome') == {str(receipt.parent / 'index')}
    assert_managed_gene_metadata(warm, receipt)
    subprocess.run(shell_command(warm, 'stage_project_annotation.sh'), cwd=work, check=True)
    local_snapshot = (work / receipt.parent / 'annotation').resolve()
    assert local_snapshot.name == annotation_snapshot(source)[0].name
    for name in ('fasta/genome.fa', 'anno/genes.gtf'):
        assert (local_snapshot / name).read_bytes() == (source / name).read_bytes()

    # A missing generation or manifest still fails while count selects that pair.
    generation = (work / receipt.parent / 'index').resolve(strict=True)
    (generation / removed).rename(tmp_path / 'withdrawn')
    damaged = plan_existing(work, COUNT_TARGET, expect_success=False)
    assert 'cellranger count ' not in damaged and str(generation) in damaged

    # Retain completed metadata from both the previous generation and the shared provider directory.
    for root in (local_snapshot, source):
        for name in ('fasta/genome.fa.fai', 'anno/exons.gtf', 'anno/exons.bed', 'fasta/exons.fa', 'anno/genes.tsv'):
            (root / name).write_text('older annotation')
    replacement = source / changed
    old = replacement.stat()
    replacement.write_bytes(replacement.read_bytes() + b'\n')
    os.utime(replacement, ns=(old.st_atime_ns, old.st_mtime_ns))

    fresh = tmp_path / 'fresh'
    fresh.mkdir()
    for plan in (plan_existing(work, 'multiqc_report'), dry_run(fresh, organism, ext_dir=str(shared), **provider_config)):
        digest = annotation_snapshot(source)[0].name
        assert arguments(plan, 'cellranger mkref ', '--expected-source') == {digest}
        new_receipt = Path(next(iter(arguments(plan, 'cellranger mkref ', '--receipt'))))
        assert new_receipt != receipt
        assert all(digest in Path(value).parts for value in arguments(plan, 'cellranger count ', '--transcriptome'))
        assert_managed_gene_metadata(plan, new_receipt)
        for rule in ('exon_gtf', 'exon_bed', 'exon_fasta', 'convert_gtf2gene_info', 'multiqc_report'):
            assert f'rule {rule}:' in plan
    assert (local_snapshot / changed).read_bytes() != replacement.read_bytes()


def test_vendor_default_release_through_reporting(tmp_path):
    log = dry_run(tmp_path, 'homo_sapiens')
    release = DB_DEFAULTS['10xgenomics']['release']
    ref = f'data/ext/10xgenomics/release-{release}/homo_sapiens/GRCh38'
    assert 'rule multiqc_report:' in log
    assert f'refdata-gex-GRCh38-{release}.tar.gz' in log
    assert arguments(log, 'cellranger count ', '--transcriptome') == {ref}
    assert_gene_metadata(log, ref)
    assert 'cellranger mkref ' not in log


@pytest.mark.parametrize('release', [None, ''])
def test_vendor_without_effective_release(tmp_path, release):
    log = dry_run(tmp_path, 'homo_sapiens', COUNT_TARGET, expect_success=False, db={'10xgenomics': {'release': release}})
    assert 'homo_sapiens' in log and 'db.10xgenomics.release' in log
    assert 'cellranger mkref ' not in log


def write_panel(work, ref, missing=None):
    for file, contents in [('fasta/genome.fa', '>synthetic\nACGT\n'), ('anno/genes.gtf', '# synthetic annotation\n')]:
        if file != missing:
            path = work / ref / file
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(contents)


@pytest.mark.parametrize('ref,organism', [('refs/panel-v1', 'synthetic_organism'), ('refs/10xgenomics/panel two', 'homo_sapiens')])
def test_project_directory_through_reporting(tmp_path, ref, organism):
    write_panel(tmp_path, ref)
    db = {'10xgenomics': {'release': None}}
    if organism == 'synthetic_organism':
        db.update(ensembl={'release': None}, ensemblgenomes={'release': None})
    log = dry_run(tmp_path, organism, ref_dir=ref, db=db)
    assert 'rule multiqc_report:' in log
    assert_build_and_count(log, ref)
    assert_gene_metadata(log, annotation_snapshot(tmp_path / ref)[0])
    assert 'refdata-gex-' not in log


def annotation_snapshot(panel):
    checksums = ''.join(f'{hashlib.sha256((panel / path).read_bytes()).hexdigest()}  {path}\n'
                        for path in ('fasta/genome.fa', 'anno/genes.gtf'))
    digest = hashlib.sha256(checksums.encode()).hexdigest()
    return Path('data/tmp/singlecell/quant/reference') / digest, checksums


@pytest.mark.parametrize('changed', ['fasta/genome.fa', 'anno/genes.gtf'])
def test_report_refreshes_gene_table_and_exons_for_changed_pair(tmp_path, changed):
    panel = tmp_path / 'panel'
    write_panel(tmp_path, panel)
    old_root, old_checksums = annotation_snapshot(panel)
    old_plan = dry_run(tmp_path, 'homo_sapiens', ref_dir=str(panel))
    assert_gene_metadata(old_plan, old_root)

    # Stand in for the previous generation and for timestamp-current legacy metadata.
    # Exercise only the staging/checksum helper on tiny text files.
    stage = ['bash', WORKFLOW.parent / 'rules/quant/scripts/stage_project_annotation.sh', panel / 'fasta/genome.fa', panel / 'anno/genes.gtf']
    subprocess.run([*stage, tmp_path / old_root, old_checksums], check=True)
    for root in (tmp_path / old_root, panel):
        for path in ('fasta/genome.fa.fai', 'anno/exons.gtf', 'anno/exons.bed', 'fasta/exons.fa', 'anno/genes.tsv'):
            (root / path).write_text('previous generation')
    source = panel / changed
    before = source.stat()
    source.write_bytes(source.read_bytes() + b'\n')
    os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))

    new_root, new_checksums = annotation_snapshot(panel)
    assert new_root != old_root
    plan = plan_existing(tmp_path, 'multiqc_report')
    for rule in ('cellranger_mkref', 'project_annotation_inputs', 'exon_gtf', 'exon_bed', 'exon_fasta', 'convert_gtf2gene_info', 'multiqc_report'):
        assert f'rule {rule}:' in plan
    assert 'rule export_project_gene_table:' not in plan
    assert all(new_root.name in Path(value).parts for value in arguments(plan, 'cellranger count ', '--transcriptome'))
    assert_gene_metadata(plan, new_root)
    assert str(old_root) not in plan

    stale = subprocess.run([*stage, tmp_path / old_root, old_checksums], capture_output=True, text=True)
    assert stale.returncode and 'FAILED' in stale.stdout
    assert (tmp_path / old_root / changed).read_bytes() != source.read_bytes()
    subprocess.run([*stage, tmp_path / new_root, new_checksums], check=True)
    assert (tmp_path / new_root / changed).read_bytes() == source.read_bytes()
    assert not (tmp_path / new_root / 'fasta/genome.fa.fai').exists()


def test_existing_annotation_snapshot_used_as_project_reference(tmp_path):
    initial = tmp_path / 'initial'
    write_panel(tmp_path, initial)
    snapshot, _ = annotation_snapshot(initial)
    selected = tmp_path / snapshot
    selected.parent.mkdir(parents=True)
    initial.rename(selected)
    (tmp_path / 'selected-panel').symlink_to(selected, target_is_directory=True)
    plan = dry_run(tmp_path, 'homo_sapiens', ref_dir='selected-panel')
    assert_gene_metadata(plan, snapshot / 'derived')
    assert 'rule export_project_gene_table:' not in plan


@pytest.mark.parametrize('workflow_file,kit,rule', [
    (WORKFLOW, TENX_KIT, 'singlecell_fastq_screen'),
    (WORKFLOW.parents[1] / 'rnaseq/rnaseq.smk', 'Lexogen SENSE Total RNA-Seq Library Prep Kit (w/RiboCop rRNA Depletion Kit V1.2)', 'fastq_screen'),
])
def test_project_genome_in_fastq_screen(tmp_path, workflow_file, kit, rule):
    ref = Path('refs/panel two')
    write_panel(tmp_path, ref)
    target = f'data/tmp/{workflow_file.stem}/qc/fastq_screen/libA_screen.txt'
    log = dry_run(tmp_path, 'synthetic_organism', target, workflow_file=workflow_file, kit=kit, ref_dir=str(ref),
                  db={'reference_db': 'ensembl', 'ensembl': {'release': None}})
    assert f'rule {rule}:' in log
    fasta = ref / 'fasta/genome.fa'
    digest = hashlib.sha256((tmp_path / fasta).read_bytes()).hexdigest()
    index_dir = Path(f'data/tmp/{workflow_file.stem}/qc/fastq_screen/indexes/{digest}')
    build = shell_command(log, 'build_project_screen_index.sh')
    assert build[2:] == [str(fasta), digest, 'bt2', '1', str(index_dir)]

    records = next(shlex.split(line) for line in log.splitlines() if 'DATABASE Synthetic_Organism ' in line)
    project_prefix = next(record.split()[2] for record in records if record.startswith('DATABASE Synthetic_Organism '))
    human_dir = tmp_path / 'data/ext/ensembl' / f"release-{DB_DEFAULTS['ensembl']['release']}" / 'homo_sapiens/GRCh38'
    assert f'DATABASE Human {human_dir}/index/genome/bowtie2/genome' in records
    assert 'release-None' not in log
    # The same per-job shell invocation installs the directory alias and starts FastQ Screen.
    screen = next(shlex.split(line) for line in log.splitlines() if '&& fastq_screen ' in line)
    pos = screen.index('ln')
    assert screen[pos:pos + 7] == ['ln', '-s', '--', str(tmp_path / index_dir / 'genome.*.bt2'), str(Path(project_prefix).parent) + '/', '&&', 'fastq_screen']
    assert Path(project_prefix).name == 'genome'


@pytest.mark.parametrize('missing', ['fasta/genome.fa', 'anno/genes.gtf'])
def test_missing_project_reference_input(tmp_path, missing):
    ref = 'refs/panel-v1'
    write_panel(tmp_path, ref, missing=missing)
    log = dry_run(tmp_path, 'homo_sapiens', COUNT_TARGET, ref_dir=ref, expect_success=False)
    assert f'{ref}/{missing}' in log
    assert 'cellranger mkref ' not in log


@pytest.mark.parametrize('provider', ['custom_reference', 'ensembl'])
def test_configured_custom_assembly(tmp_path, provider):
    db = {'reference_db': provider, 'custom_reference': {'assembly': 'Rtarandus_2.1'}}
    log = dry_run(tmp_path, 'rangifer_tarandus', COUNT_TARGET, db=db)
    assert_build_and_count(log, 'data/ext/custom_reference/rangifer_tarandus/Rtarandus_2.1')
    if provider == 'ensembl':
        assert re.search(r'no hits in ensembl lookup file.*custom_reference', log, re.I)


@pytest.mark.parametrize('changed', ['fasta/genome.fa', 'anno/genes.gtf', 'assembly', 'image'])
def test_retained_project_index_requires_current_build_identity(tmp_path, changed):
    ref = tmp_path / 'panel'
    write_panel(tmp_path, ref)
    # A complete legacy index, newer than both inputs, cannot establish which inputs and tool built it.
    for name in ('reference.json', 'fasta/genome.fa', 'genes/genes.gtf.gz', 'anno/genes.gtf'):
        path = ref / 'cellranger' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('stock')

    def planned_index(project, assembly='Panel.v1', image='example/cellranger:1'):
        work = tmp_path / project
        work.mkdir()
        log = dry_run(work, 'homo_sapiens', COUNT_TARGET, ref_dir=str(ref),
                      db={'10xgenomics': {'assembly': assembly}}, docker={'cellranger': image})
        assert_build_and_count(log, str(ref))
        directory, = arguments(log, 'cellranger count ', '--transcriptome')
        return Path(directory)

    original = planned_index('legacy')
    record = ''.join(f'{hashlib.sha256((ref / path).read_bytes()).hexdigest()}  {path}\n'
                     for path in ('fasta/genome.fa', 'anno/genes.gtf'))
    assert original.parent == ref / 'cellranger/builds' / hashlib.sha256(record.encode()).hexdigest()
    # Same normalized Cell Ranger name and image can share a generation across fresh workdirs.
    assert planned_index('same', assembly='Panel_v1') == original

    settings = {}
    if changed == 'assembly':
        settings['assembly'] = 'Panel.v2'
    elif changed == 'image':
        settings['image'] = 'example/cellranger:2'
    else:
        source = ref / changed
        old = source.stat()
        source.write_bytes(source.read_bytes() + b'\n')
        os.utime(source, ns=(old.st_atime_ns, old.st_mtime_ns))
    replacement = planned_index('changed', **settings)
    assert replacement != original
    assert (replacement.parent == original.parent) == (changed in ('assembly', 'image'))


def test_retained_receipt_checks_index_before_count(tmp_path):
    panel = tmp_path / 'panel'
    write_panel(tmp_path, panel)
    cold = dry_run(tmp_path, 'homo_sapiens', COUNT_TARGET, ref_dir=str(panel))
    build = next(shlex.split(line.removeprefix('Shell command: ')) for line in cold.splitlines()
                 if line.startswith('Shell command: ') and 'cellranger mkref ' in line)
    # Issue a real receipt using the publisher and synthetic mkref output.
    env = protected_env(tmp_path)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    subprocess.run([sys.executable, *build[1:]], cwd=tmp_path, env=env, check=True)

    warm = plan_existing(tmp_path, COUNT_TARGET)
    assert 'cellranger mkref ' not in warm
    assert 'cellranger count ' in warm
    receipt = tmp_path / build[build.index('--receipt') + 1]
    manifest = receipt.resolve(strict=True)

    # Switching to an identically prepared directory must honor the newly selected cache.
    other = tmp_path / 'second panel'
    write_panel(tmp_path, other)
    configfile = tmp_path / 'config.yaml'
    original_config = configfile.read_text()
    config = yaml.safe_load(original_config)
    config['ref_dir'] = str(other)
    configfile.write_text(yaml.safe_dump(config))
    switched = plan_existing(tmp_path, COUNT_TARGET)
    assert_build_and_count(switched, str(other))
    other_receipt, = arguments(switched, 'cellranger mkref ', '--receipt')
    assert (tmp_path / other_receipt) != receipt
    assert receipt.resolve(strict=True) == manifest

    alias = tmp_path / 'first panel alias'
    alias.symlink_to(panel, target_is_directory=True)
    config['ref_dir'] = str(alias)
    configfile.write_text(yaml.safe_dump(config))
    same = plan_existing(tmp_path, COUNT_TARGET)
    assert 'cellranger mkref ' not in same
    transcriptome, = arguments(same, 'cellranger count ', '--transcriptome')
    assert Path(transcriptome).resolve() == manifest.parent
    configfile.write_text(original_config)

    for removed in ('fasta/genome.fa', 'star/genomeParameters.txt', 'star/sjdbInfo.txt'):
        member = manifest.parent / removed
        original = member.read_bytes()
        member.unlink()
        damaged = plan_existing(tmp_path, COUNT_TARGET, expect_success=False)
        assert 'cellranger count ' not in damaged and str(member) in damaged
        member.write_bytes(original)

    annotation = receipt.parent / 'annotation'
    target = annotation.readlink()
    annotation.unlink()
    annotation.symlink_to(tmp_path)  # An existing directory for a different annotation generation.
    damaged = plan_existing(tmp_path, COUNT_TARGET, expect_success=False)
    assert 'cellranger count ' not in damaged and 'Invalid Cell Ranger annotation pin' in damaged

    annotation.unlink()
    annotation.symlink_to(target)
    # The same receipt must also agree with the generation's recorded identity.
    manifest.write_text('{"genome": "another_panel"}')
    damaged = plan_existing(tmp_path, COUNT_TARGET, expect_success=False)
    assert 'cellranger count ' not in damaged and str(manifest) in damaged
