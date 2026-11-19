"""Plan the BFQ reporting target and inspect the tools' reference arguments with an empty cache."""
import re
import shlex
import subprocess
import sys

import pytest
import yaml

from test_velocity import TENX_KIT, WORKFLOW, protected_env

DB_DEFAULTS = yaml.safe_load((WORKFLOW.parents[1] / 'main.config').read_text())['db']
COUNT_TARGET = '_tmp_libA/outs/web_summary.html'


def dry_run(work, organism, target='multiqc_report', expect_success=True, **settings):
    config = dict(workflow='singlecell', organism=organism, project_id=['SYNTHETIC'], libprepkit=TENX_KIT,
                  machine='Synthetic sequencer', read_geometry=[28, 90],
                  samples={'libA': {'R1': 'libA_R1.fastq.gz', 'R2': 'libA_R2.fastq.gz'}}, **settings)
    (work / 'config.yaml').write_text(yaml.safe_dump(config))
    (work / 'pep.yaml').write_text('pep_version: 2.0.0\nsample_table: samples.csv\n')
    (work / 'samples.csv').write_text('sample_name,R1,R2,paired_end\nlibA,libA_R1.fastq.gz,libA_R2.fastq.gz,1\n')
    fastqs = work / 'data/raw/fastq'
    fastqs.mkdir(parents=True)
    for read in config['samples']['libA'].values():
        (fastqs / read).touch()
    (work / 'Snakefile').write_text(f'pepfile: "pep.yaml"\nconfigfile: "config.yaml"\ninclude: {str(WORKFLOW)!r}\n')
    assert not (work / 'data/ext').exists()
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
    assert f'{ref}/cellranger/reference.json' in log
    assert arguments(log, 'cellranger count ', '--transcriptome') == {f'{ref}/cellranger'}


def assert_gene_metadata(log, ref):
    assert arguments(log, 'bedtools getfasta ', '-fi') == {f'{ref}/fasta/genome.fa'}
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
    assert_gene_metadata(log, ref)
    assert 'rule txgenomics_org_prebuild:' not in log
    assert 'refdata-gex-' not in log


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
    assert_gene_metadata(log, ref)
    assert 'refdata-gex-' not in log


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
