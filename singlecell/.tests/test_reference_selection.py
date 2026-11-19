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


def plan_existing(work, target, expect_success=True, extra=()):
    result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', *extra, target],
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


def assert_build_and_count(log, ref, work):
    assert arguments(log, 'cellranger mkref ', '--fasta') == {f'{ref}/fasta/genome.fa'}
    assert arguments(log, 'cellranger mkref ', '--genes') == {f'{ref}/anno/genes.gtf'}
    assert arguments(log, 'cellranger mkref ', '--cache') == {str((work / ref / 'cellranger/builds').resolve())}
    receipt, = arguments(log, 'cellranger mkref ', '--receipt')
    assert arguments(log, 'cellranger count ', '--transcriptome') == {str(Path(receipt).parent / 'index')}


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
    assert_build_and_count(log, ref, tmp_path)
    # Keep the assembly in the paths and provenance, but pass a valid Cell Ranger genome identifier.
    assert arguments(log, 'cellranger mkref ', '--genome') == {assembly.replace('.', '_')}
    receipt, = arguments(log, 'cellranger mkref ', '--receipt')
    assert_managed_gene_metadata(log, Path(receipt))
    assert 'rule txgenomics_org_prebuild:' not in log
    assert 'refdata-gex-' not in log


def test_ensembl_genomes_sources_for_direct_and_10x_routes(tmp_path):
    organism, assembly = 'arabidopsis_thaliana', 'TAIR10'
    release = DB_DEFAULTS['ensemblgenomes']['release']
    old_release = str(int(release) - 1)
    main_release = DB_DEFAULTS['ensembl']['release']
    shared = tmp_path / 'shared'
    # The legacy directory has both files, but its path does not distinguish Plants releases.
    legacy = shared / 'ensembl' / f'release-{main_release}' / organism / assembly
    write_panel(tmp_path, legacy)
    work = tmp_path / 'retained'
    work.mkdir()

    def check_cold(log, directory, selected):
        archive = f'https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-{selected}/plants'
        genome_url = f'{archive}/fasta/{organism}/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'
        gtf_url = f'{archive}/gtf/{organism}/Arabidopsis_thaliana.TAIR10.{selected}.gtf.gz'
        assert f'-O- {genome_url} |' in log and f'-O- {gtf_url} |' in log
        assert f'Ensembl DNA,release-{selected},{genome_url}' in log and f'Ensembl GTF,release-{selected},{gtf_url}' in log
        location = Path('ensemblgenomes') / f'release-{selected}' / 'plants' / organism / assembly
        ref = shared / location
        assert str(ref / 'anno' / f'{organism}.{assembly}.{selected}.gtf') in log
        assert_build_and_count(log, ref, directory)
        receipt, = arguments(log, 'cellranger mkref ', '--receipt')
        assert f'/{location}/' in receipt
        assert_managed_gene_metadata(log, Path(receipt))
        assert 'rule multiqc_report:' in log and 'refdata-gex-' not in log
        return ref

    # A cold fallback fetches the selected Plants release. Direct Ensembl can reuse that same pair.
    fallback = dry_run(work, organism, ext_dir=str(shared))
    assert f'source=ensemblgenomes release={release} assembly={assembly}' in fallback and 'division=plants' in fallback
    ref = check_cold(fallback, work, release)
    write_panel(tmp_path, ref)
    configfile = work / 'config.yaml'
    config = yaml.safe_load(configfile.read_text())
    config['db'] = {'reference_db': 'ensembl'}
    configfile.write_text(yaml.safe_dump(config))
    warm = plan_existing(work, 'multiqc_report')
    assert 'rule ensembl_genome:' not in warm and 'rule ensembl_gtf:' not in warm
    assert_build_and_count(warm, ref, work)

    # Changing only the division release must retrieve both archives on a retained or a fresh run.
    config['db']['ensemblgenomes'] = {'release': old_release}
    configfile.write_text(yaml.safe_dump(config))
    check_cold(plan_existing(work, 'multiqc_report'), work, old_release)
    fresh = tmp_path / 'fresh'
    fresh.mkdir()
    check_cold(dry_run(fresh, organism, ext_dir=str(shared), db=config['db']), fresh, old_release)


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


def retain_count_outputs(work, plan, env):
    # Retain completed count outputs and their real Snakemake provenance without running Cell Ranger.
    outputs = re.search(r'^rule cellranger_quant_:\n    input: .+\n    output: (.+)$', plan, re.MULTILINE).group(1)
    for name in outputs.split(', '):
        path = work / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    subprocess.run([sys.executable, '-m', 'snakemake', '--cores', '1', '--touch', '--force', COUNT_TARGET],
                   cwd=work, env=env, check=True, capture_output=True, text=True)
    assert 'Nothing to be done' in plan_existing(work, COUNT_TARGET)


def test_forced_warm_publisher_keeps_previously_planned_source(tmp_path):
    cold = dry_run(tmp_path, 'salmo_salar', COUNT_TARGET)
    initial = shell_command(cold, 'cellranger mkref ')
    receipt = Path(initial[initial.index('--receipt') + 1])
    assert '--expected-source' not in initial
    # Have Snakemake record the real cold publisher's inputs and params. Only the tiny
    # synthetic getter and the publisher using the mock mkref executable may run.
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text(snakefile.read_text() + """
ruleorder: synthetic_provider > ensembl_genome
ruleorder: synthetic_provider > ensembl_symlink_genes
rule synthetic_provider:
    output: get_genome(), get_gtf()
    shell: "echo '>synthetic' > {output[0]:q}; echo 'ACGT' >> {output[0]:q}; echo '# synthetic' > {output[1]:q}"
""")
    env = protected_env(tmp_path)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    cmd = [sys.executable, '-m', 'snakemake', '--cores', '1', '-s', str(snakefile), str(receipt),
           '--allowed-rules', 'synthetic_provider', 'cellranger_mkref']
    subprocess.run(cmd, cwd=tmp_path, env=env, check=True, capture_output=True, text=True)

    warm = plan_existing(tmp_path, COUNT_TARGET)
    assert 'rule cellranger_mkref:' not in warm
    retain_count_outputs(tmp_path, warm, env)
    forced = plan_existing(tmp_path, str(receipt), extra=['--force'])
    publish = shell_command(forced, 'cellranger mkref ')
    local = tmp_path / receipt.parent
    generation = (local / 'index').resolve(strict=True)
    manifest = (local / 'ready').read_bytes()
    annotation = (local / 'annotation').resolve(strict=True)
    gtf = tmp_path / publish[publish.index('--genes') + 1]
    original = gtf.read_bytes()
    gtf.write_bytes(original + b'# changed after planning\n')
    failed = subprocess.run(publish, cwd=tmp_path, env=env, capture_output=True, text=True)
    assert failed.returncode and 'changed after planning' in failed.stderr, failed.stderr
    assert publish[publish.index('--expected-source') + 1] == generation.parent.name
    assert (local / 'index').resolve(strict=True) == generation
    assert (local / 'ready').read_bytes() == manifest
    assert (local / 'annotation').resolve(strict=True) == annotation
    assert (generation / 'gcf-source/anno/genes.gtf').read_bytes() == original


def test_managed_refetch_after_partial_cache_loss(tmp_path):
    work = tmp_path / 'work'
    work.mkdir()
    ext = tmp_path / 'ext'
    cold = dry_run(work, 'salmo_salar', ext_dir=str(ext))
    publish = shell_command(cold, 'cellranger mkref ')
    source = Path(publish[publish.index('--source-dir') + 1])
    write_panel(tmp_path, source)
    env = protected_env(work)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    subprocess.run([sys.executable, *publish[1:]], cwd=work, env=env, check=True)
    receipt = Path(publish[publish.index('--receipt') + 1])
    fasta = source / 'fasta/genome.fa'
    gtf = source / 'anno/genes.gtf'

    for run in range(2):
        # The first obsolete generation is gone. The next is intact, and the refetch
        # returns the same bytes. Each pending download must precede reference selection.
        generation = (work / receipt.parent / 'index').resolve(strict=True)
        if run == 0:
            generation.rename(tmp_path / 'old-index')
            old = gtf.stat()
            gtf.write_bytes(gtf.read_bytes() + b'\n')
            os.utime(gtf, ns=(old.st_atime_ns, old.st_mtime_ns))
        content = fasta.read_bytes()
        fasta.unlink()
        recovery = plan_existing(work, 'multiqc_report')
        assert 'rule ensembl_genome:' in recovery
        publish = shell_command(recovery, 'cellranger mkref ')
        new_receipt = Path(publish[publish.index('--receipt') + 1])
        assert new_receipt != receipt
        assert '--expected-source' not in publish
        assert arguments(recovery, 'cellranger count ', '--transcriptome') == {str(new_receipt.parent / 'index')}
        assert_managed_gene_metadata(recovery, new_receipt)
        assert 'rule multiqc_report:' in recovery

        # Stand in for the getter and mkref, retaining the provider's old mtimes.
        fasta.write_bytes(content)
        os.utime(fasta, (100, 100))
        subprocess.run([sys.executable, *publish[1:]], cwd=work, env=env, check=True)
        warm = plan_existing(work, 'multiqc_report')
        assert 'rule cellranger_mkref:' not in warm
        assert arguments(warm, 'cellranger count ', '--transcriptome') == {str(new_receipt.parent / 'index')}
        assert_managed_gene_metadata(warm, new_receipt)
        subprocess.run(shell_command(warm, 'stage_project_annotation.sh'), cwd=work, check=True)
        local = work / new_receipt.parent / 'annotation'
        assert (local / 'anno/genes.gtf').read_bytes() == gtf.read_bytes()
        for name in ('fasta/genome.fa.fai', 'anno/exons.gtf', 'anno/exons.bed', 'fasta/exons.fa', 'anno/genes.tsv'):
            (local / name).write_text('previous result')
        if run == 0:
            retain_count_outputs(work, warm, env)
        receipt = new_receipt

    forced = plan_existing(work, str(receipt), extra=['--force'])
    assert arguments(forced, 'cellranger mkref ', '--expected-source') == {annotation_snapshot(source)[0].name}


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
    assert_build_and_count(log, ref, tmp_path)
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
    stage = ['bash', WORKFLOW.parents[1] / 'gcfdb/scripts/stage_project_annotation.sh', panel / 'fasta/genome.fa', panel / 'anno/genes.gtf']
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
    stock_name = 'Homo_sapiens.GRCh38.dna.primary_assembly'
    human_dir = tmp_path / 'data/ext/fastq_screen/ensembl' / f"release-{DB_DEFAULTS['ensembl']['release']}" / 'homo_sapiens/GRCh38'
    assert f'DATABASE Human {human_dir}/index/{stock_name}/bowtie2/{stock_name}' in records
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
    assert_build_and_count(log, 'data/ext/custom_reference/rangifer_tarandus/Rtarandus_2.1', tmp_path)
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
        assert_build_and_count(log, ref, work)
        command = shell_command(log, 'cellranger mkref ')
        env = protected_env(work)
        env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
        subprocess.run([sys.executable, *command[1:]], cwd=work, env=env, check=True)
        transcriptome, = arguments(log, 'cellranger count ', '--transcriptome')
        return (work / transcriptome).resolve(strict=True)

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


def test_selected_project_symlink_stays_pinned_after_planning(tmp_path):
    panel, other, alias = (tmp_path / name for name in ('panel', 'other', 'current'))
    write_panel(tmp_path, panel)
    write_panel(tmp_path, other)
    alias.symlink_to(panel, target_is_directory=True)
    log = dry_run(tmp_path, 'homo_sapiens', ref_dir=str(alias))
    assert_build_and_count(log, panel, tmp_path)
    command = shell_command(log, 'cellranger mkref ')
    alias.unlink()
    alias.symlink_to(other, target_is_directory=True)
    env = protected_env(tmp_path)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    subprocess.run([sys.executable, *command[1:]], cwd=tmp_path, env=env, check=True)
    transcriptome, = arguments(log, 'cellranger count ', '--transcriptome')
    generation = (tmp_path / transcriptome).resolve(strict=True)
    assert generation.is_relative_to(panel / 'cellranger/builds')
    assert not (other / 'cellranger').exists()
    metadata = next(line for line in log.splitlines() if line.startswith('Shell command:') and 'stage_project_annotation.sh' in line)
    assert shlex.join([str(panel / 'fasta/genome.fa'), str(panel / 'anno/genes.gtf')]) in metadata


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
    assert_build_and_count(switched, other, tmp_path)
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
    alias.unlink()
    alias.symlink_to(other, target_is_directory=True)
    assert (tmp_path / transcriptome).resolve(strict=True) == manifest.parent
    configfile.write_text(original_config)

    for removed in ('fasta/genome.fa', 'star/genomeParameters.txt', 'star/sjdbInfo.txt'):
        member = manifest.parent / removed
        original = member.read_bytes()
        member.unlink()
        damaged = plan_existing(tmp_path, COUNT_TARGET, expect_success=False)
        assert 'cellranger count ' not in damaged and str(member) in damaged
        member.write_bytes(original)

    for name, truncated in (('reference.json', True), ('star/SA', False)):
        member = manifest.parent / name
        original = member.read_bytes()
        old = member.stat()
        member.write_bytes(b'' if truncated else b'x' * len(original))
        os.utime(member, ns=(old.st_atime_ns, old.st_mtime_ns))
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


def test_only_jobs_using_cellranger_pin_validate_it(tmp_path):
    cold = dry_run(tmp_path, 'salmo_salar')
    build = shell_command(cold, 'cellranger mkref ')
    write_panel(tmp_path, Path(build[build.index('--source-dir') + 1]))
    receipt = Path(build[build.index('--receipt') + 1])
    env = protected_env(tmp_path)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    subprocess.run([sys.executable, *build[1:]], cwd=tmp_path, env=env, check=True)

    # Count and managed annotation share one full verification in the report DAG.
    with (tmp_path / 'Snakefile').open('a') as rules:
        rules.write("""
verify_original = verify_cellranger_receipt
def verify_cellranger_receipt(*args):
    with open('receipt-checks', 'a') as checks:
        print(args[0], file=checks)
    return verify_original(*args)
""")
    checks = tmp_path / 'receipt-checks'
    assert 'rule multiqc_report:' in plan_existing(tmp_path, 'multiqc_report')
    assert checks.read_text().splitlines() == [str(receipt)]

    (tmp_path / receipt.parent / 'index/star/SA').unlink()
    annotations = str(receipt.parent / 'annotation/anno/genes.tsv')
    for target in [COUNT_TARGET, annotations, 'multiqc_report']:
        damaged = plan_existing(tmp_path, target, expect_success=False)
        assert 'Incomplete Cell Ranger reference' in damaged and '/star/SA' in damaged
    checks.unlink()
    qc = plan_existing(tmp_path, 'data/tmp/singlecell/qc/fastqc/libA.fastqc.zip')
    assert 'fastqc -t' in qc and not checks.exists()

    settings = tmp_path / 'config.yaml'
    config = yaml.safe_load(settings.read_text())
    config['quant'] = {'method': '10x_starsolo', 'aggregate': {'method': 'scanpy'}, 'use_velo': False}
    settings.write_text(yaml.safe_dump(config))
    star = plan_existing(tmp_path, 'multiqc_report')
    assert 'rule starsolo_quant:' in star and 'rule multiqc_report:' in star
    assert 'cellranger count ' not in star and not checks.exists()


def test_mixed_methods_select_the_counting_reference_for_annotation(tmp_path):
    quant = {'method': 'cellranger,10x_starsolo', 'aggregate': {'method': 'scanpy'}, 'use_velo': False}
    source = Path(f"data/ext/ensembl/release-{DB_DEFAULTS['ensembl']['release']}/salmo_salar/Ssal_v3.1")
    provider_table = str(source / 'anno/genes.tsv')
    cold = dry_run(tmp_path, 'salmo_salar', '_tmp/10x_starsolo/raw/libA/anndata.light.h5ad', quant=quant)
    assert 'rule starsolo_quant:' in cold and 'rule cellranger_mkref:' not in cold
    assert_gene_metadata(cold, source)

    combined = plan_existing(tmp_path, 'quant_all')
    build = shell_command(combined, 'cellranger mkref ')
    receipt = Path(build[build.index('--receipt') + 1])
    assert arguments(combined, 'cellranger count ', '--transcriptome') == {str(receipt.parent / 'index')}
    tables = set()
    for line in combined.splitlines():
        if line.startswith('Shell command:') and 'convert_scanpy.py ' in line:
            args = shlex.split(line)
            tables.add((args[args.index('-f') + 1], args[args.index('--feature-info') + 1]))
    assert tables == {('cellranger', str(receipt.parent / 'annotation/anno/genes.tsv')), ('10x_starsolo', provider_table)}

    write_panel(tmp_path, source)
    env = protected_env(tmp_path)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    subprocess.run([sys.executable, *build[1:]], cwd=tmp_path, env=env, check=True)
    (tmp_path / receipt.parent / 'index/star/SA').unlink()
    filtered = plan_existing(tmp_path, '_tmp/10x_starsolo/filtered/libA/anndata.light.h5ad')
    assert 'rule starsolo_quant:' in filtered and 'rule cellranger_mkref:' not in filtered
    assert_gene_metadata(filtered, source)
    damaged = plan_existing(tmp_path, '_tmp/cellranger/filtered/libA/anndata.light.h5ad', expect_success=False)
    assert 'Incomplete Cell Ranger reference' in damaged and '/star/SA' in damaged


def test_mixed_methods_share_explicit_reference_metadata(tmp_path):
    panel = Path('refs/augmented panel')
    write_panel(tmp_path, panel)
    quant = {'method': 'cellranger,10x_starsolo', 'aggregate': {'method': 'scanpy'}, 'use_velo': False}
    plan = dry_run(tmp_path, 'salmo_salar', 'quant_all', ref_dir=str(panel), quant=quant)
    assert 'rule starsolo_quant:' in plan and 'cellranger count ' in plan
    assert_gene_metadata(plan, annotation_snapshot(tmp_path / panel)[0])


@pytest.mark.parametrize('methods,provider,index_type,changed', [
    ('cellranger,10x_starsolo', '10xgenomics', 'starsolo', 'anno/genes.gtf'),
    ('10x_starsolo', 'ensembl', 'star', 'fasta/genome.fa'),
])
def test_starsolo_recounts_after_project_reference_changes(tmp_path, methods, provider, index_type, changed):
    panel = Path('refs/panel-v1')
    write_panel(tmp_path, panel)
    target = '_tmp/10x_starsolo/filtered/libA/anndata.light.h5ad'
    cold = dry_run(tmp_path, 'salmo_salar', target, ref_dir=str(panel), db={'reference_db': provider},
                   quant={'method': methods, 'use_velo': False})
    env = protected_env(tmp_path)

    def stage_inputs(snapshot):
        # Only the workflow's checksum-verifying copy rule runs on these synthetic sequences.
        subprocess.run([sys.executable, '-m', 'snakemake', '--cores', '1', str(snapshot / 'fasta/genome.fa'),
                        '--allowed-rules', 'project_annotation_inputs'], cwd=tmp_path, env=env, check=True, capture_output=True, text=True)
        for member in ('fasta/genome.fa', 'anno/genes.gtf'):
            assert (tmp_path / snapshot / member).read_bytes() == (tmp_path / panel / member).read_bytes()

    original, _ = annotation_snapshot(tmp_path / panel)
    stage_inputs(original)
    # Give the complete synthetic DAG real Snakemake provenance without executing any scientific command.
    for outputs in re.findall(r'^    output: (.+)$', cold, re.MULTILINE):
        for name in outputs.split(', '):
            path = tmp_path / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
    subprocess.run([sys.executable, '-m', 'snakemake', '--cores', '1', '--touch', '--forceall', '--notemp', target],
                   cwd=tmp_path, env=env, check=True, capture_output=True, text=True)
    assert 'Nothing to be done' in plan_existing(tmp_path, target)

    path = tmp_path / panel / changed
    previous = path.read_bytes()
    stat = path.stat()
    path.write_bytes(previous + b'\n')
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    current, _ = annotation_snapshot(tmp_path / panel)
    retained = plan_existing(tmp_path, target)
    assert current != original
    assert 'rule starsolo_quant:' in retained and 'rule tmp_lightweight_filtered:' in retained
    for plan, snapshot in ((cold, original), (retained, current)):
        index_dir = snapshot / 'index/genome' / index_type / 'r_90'
        assert arguments(plan, '--runMode genomeGenerate', '--genomeFastaFiles') == {str(snapshot / 'fasta/genome.fa')}
        assert arguments(plan, '--runMode genomeGenerate', '--sjdbGTFfile') == {str(snapshot / 'anno/genes.gtf')}
        assert arguments(plan, '--runMode genomeGenerate', '--genomeDir') == {str(index_dir)}
        index_file = str(index_dir / 'SA')
        assert f'    output: {index_file}\n' in plan
        count_inputs = re.search(r'^rule starsolo_quant:\n    input: (.+)$', plan, re.MULTILINE).group(1)
        assert index_file in count_inputs.split(', ')
        assert arguments(plan, 'STAR --soloCBwhitelist ', '--genomeDir') == {str(index_dir)}
        assert arguments(plan, 'convert_scanpy.py ', '--feature-info') == {str(snapshot / 'anno/genes.tsv')}
        assert 'cellranger mkref ' not in plan
    stage_inputs(current)
    assert (tmp_path / original / changed).read_bytes() == previous


@pytest.mark.parametrize('use_alias,new_recipe', [(False, False), (True, True)])
def test_project_reusing_published_cellranger_source_preserves_both_pins(tmp_path, use_alias, new_recipe):
    owner, importer = tmp_path / 'owner', tmp_path / 'importer'
    owner.mkdir()
    importer.mkdir()
    panel = tmp_path / 'panel'
    write_panel(tmp_path, panel)

    def expected_metadata(source):
        snapshot, _ = annotation_snapshot(source)
        if source == panel:
            return snapshot
        source_key = hashlib.sha256(os.fsencode(source.resolve())).hexdigest()
        return snapshot.parent / 'imported' / source_key / snapshot.name

    def publish(work, source, genome, image):
        plan = dry_run(work, 'homo_sapiens', ref_dir=str(source), db={'10xgenomics': {'assembly': genome}},
                       docker={'cellranger': image})
        command = shell_command(plan, 'cellranger mkref ')
        env = protected_env(work)
        env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
        subprocess.run([sys.executable, *command[1:]], cwd=work, env=env, check=True)
        transcriptome, = arguments(plan, 'cellranger count ', '--transcriptome')
        assert_gene_metadata(plan, expected_metadata(source))
        return (work / transcriptome).resolve(strict=True)

    generation = publish(owner, panel, 'Panel.v1', 'example/cellranger:1')
    manifest = (generation / 'gcf-mkref.json').read_bytes()
    frozen = generation / 'gcf-source'
    selected = frozen
    if use_alias:
        selected = importer / 'frozen-panel'
        selected.symlink_to(frozen, target_is_directory=True)
    imported = publish(importer, selected, 'Other.v2' if new_recipe else 'Panel.v1',
                       'example/cellranger:2' if new_recipe else 'example/cellranger:1')
    if new_recipe:
        assert imported.is_relative_to(importer / 'data/tmp/singlecell/quant/reference/cellranger/builds')
    else:
        assert imported == generation
    assert not (frozen / 'cellranger').exists()

    for work, expected in ((owner, generation), (importer, imported)):
        retained = plan_existing(work, 'multiqc_report')
        assert 'rule cellranger_mkref:' not in retained
        pin, = arguments(retained, 'cellranger count ', '--transcriptome')
        assert (work / pin).resolve(strict=True) == expected
        assert 'rule multiqc_report:' in retained

    for root in {selected, frozen}:
        blocked = plan_existing(importer, str(root / 'anno/genes.tsv'), expect_success=False)
        assert 'published Cell Ranger generation' in blocked
    if not new_recipe:
        # A forced publisher replaces the local receipt symlink to the generation's manifest.
        assert 'rule cellranger_mkref:' in plan_existing(importer, COUNT_TARGET, extra=('--forceall',))
        alias = importer / 'generation-alias'
        alias.symlink_to(generation, target_is_directory=True)
        member = importer / 'member-alias'
        member.symlink_to(generation / 'star/SA')
        with (importer / 'Snakefile').open('a') as rules:
            rules.write("""
rule write_custom_product:
    output: 'generation-alias/new/product'
    shell: 'touch {output:q}'
rule write_linked_member:
    output: 'member-alias'
    shell: 'touch {output:q}'
rule read_frozen_source:
    input: 'generation-alias/gcf-source/fasta/genome.fa'
    output: 'local-copy.fa'
    shell: 'cp {input:q} {output:q}'
""")
        assert 'rule read_frozen_source:' in plan_existing(importer, 'local-copy.fa')
        for target in ('generation-alias/gcf-source/anno/exons.gtf', 'member-alias'):
            blocked = plan_existing(importer, target, expect_success=False, extra=('--force',))
            assert 'Cannot write inside a published Cell Ranger generation' in blocked and target in blocked
        # A named user rule registered after the workflow must also fail before output preparation in real execution.
        actual = subprocess.run([sys.executable, '-m', 'snakemake', '--cores', '1', 'write_custom_product'],
                                cwd=importer, env=protected_env(importer), capture_output=True, text=True, timeout=120)
        assert actual.returncode and 'Cannot write inside a published Cell Ranger generation' in actual.stdout + actual.stderr
        assert not (generation / 'new').exists() and member.is_symlink()
    filtered = plan_existing(importer, 'cellranger_gtf')
    assert f'{expected_metadata(frozen)}/anno/genes.gtf.filtered' in filtered

    assert sorted(p.relative_to(frozen).as_posix() for p in frozen.rglob('*') if p.is_file()) == ['anno/genes.gtf', 'fasta/genome.fa']
    assert (generation / 'gcf-mkref.json').read_bytes() == manifest


def test_splitpipe_imports_frozen_source_without_modifying_the_generation(tmp_path):
    from test_velocity import workflow_fixture

    owner, importer = tmp_path / 'owner', tmp_path / 'importer'
    owner.mkdir()
    importer.mkdir()
    panel = tmp_path / 'panel'
    write_panel(tmp_path, panel)
    exon = 'synthetic\tfixture\texon\t1\t4\t.\t+\t.\tgene_id "GENE"; transcript_id "TX";\n'
    (panel / 'anno/genes.gtf').write_text(exon)
    original = dry_run(owner, 'homo_sapiens', COUNT_TARGET, ref_dir=str(panel))
    env = protected_env(owner)
    env['PATH'] = str(MOCK_TOOLS) + os.pathsep + env['PATH']
    subprocess.run(shell_command(original, 'cellranger mkref '), cwd=owner, env=env, check=True)
    pin, = arguments(original, 'cellranger count ', '--transcriptome')
    generation = (owner / pin).resolve(strict=True)
    inventory = set(generation.rglob('*'))
    manifest = (generation / 'gcf-mkref.json').read_bytes()
    frozen = generation / 'gcf-source'

    workflow_fixture(importer, methods='splitpipe', use_velo=False)
    (importer / 'selected').symlink_to(frozen, target_is_directory=True)
    config = importer / 'config.yaml'
    config.write_text(yaml.safe_dump({**yaml.safe_load(config.read_text()), 'ref_dir': 'selected'}))
    report = plan_existing(importer, 'multiqc_report')
    assert 'rule multiqc_report:' in report and 'cellranger mkref ' not in report
    genome, = arguments(report, 'split-pipe --mode mkref ', '--fasta')
    source = Path(genome).parents[1]
    assert (importer / source).resolve().is_relative_to(importer / 'data/tmp/singlecell/quant/reference')
    assert source.name == annotation_snapshot(frozen)[0].name
    index = source / 'index/genome/splitpipe'
    assert genome == str(source / 'fasta/genome.fa')
    assert arguments(report, 'split-pipe --mode mkref ', '--genes') == {str(source / 'anno/genes.gtf')}
    assert arguments(report, 'split-pipe --mode mkref ', '--output_dir') == {str(index)}
    assert arguments(report, 'split-pipe --mode all ', '--genome_dir') == {str(index)}
    assert_gene_metadata(report, source)
    for paths in re.findall(r'^    (?:output|log|benchmark): (.+)$', report, re.MULTILINE):
        for path in paths.split(', '):
            assert not (importer / path).resolve().is_relative_to(generation), path

    # Exercise the staging and exon rules and the real index wrapper, substituting a tiny text file for the Split-Pipe index.
    tools = importer / 'tools'
    tools.mkdir()
    stub = tools / 'split-pipe'
    stub.write_text("""#!/usr/bin/env python3
from pathlib import Path
import sys
assert sys.argv[sys.argv.index('--mode') + 1] == 'mkref'
Path(sys.argv[sys.argv.index('--output_dir') + 1], 'SA').write_text('synthetic index')
""")
    stub.chmod(0o755)
    env = protected_env(importer)
    env['PATH'] = str(tools) + os.pathsep + env['PATH']
    # Supply the single extracted exon so the normal gene-table converter can also run without bedtools.
    exon_fasta = importer / source / 'fasta/exons.fa'
    exon_fasta.parent.mkdir(parents=True, exist_ok=True)
    exon_fasta.write_text('>synthetic:0-4\nACGT\n')
    command = [sys.executable, '-m', 'snakemake', '--cores', '1', str(index / 'SA'), str(source / 'anno/exons.gtf'),
               str(source / 'anno/genes.tsv'), '--allowed-rules', 'stage_frozen_reference', 'splitpipe_index',
               'splitpipe_mkref_parfile', 'exon_gtf', 'convert_gtf2gene_info']
    result = subprocess.run(command, cwd=importer, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    for member in ('fasta/genome.fa', 'anno/genes.gtf'):
        assert (importer / source / member).read_bytes() == (frozen / member).read_bytes()
    assert (importer / index / 'SA').read_text() == 'synthetic index'
    assert (importer / source / 'anno/exons.gtf').read_text() == exon
    assert (importer / source / 'anno/genes.tsv').read_text().splitlines()[1].split('\t')[0] == 'GENE'
    blocked = plan_existing(importer, 'selected/anno/exons.gtf', expect_success=False)
    assert 'Cannot write inside a published Cell Ranger generation' in blocked
    assert set(generation.rglob('*')) == inventory
    assert (generation / 'gcf-mkref.json').read_bytes() == manifest
    retained = plan_existing(owner, COUNT_TARGET)
    assert arguments(retained, 'cellranger count ', '--transcriptome') == {pin}
    assert 'cellranger mkref ' not in retained
