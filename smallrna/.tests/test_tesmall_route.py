"""TEsmall route (#374): selection, database planning on the mirflowz reference (one database per UNITAS dump, the
genome shared), run settings, guards and the joined table."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import pytest

from test_quant_methods import MIRBASE_GFF, MIRBASE_GFF_ID, configure, dry_run, mirbase_gff_id, rule_block, snakemake

ROOT = Path(__file__).resolve().parents[2]
REF = 'data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus'
GENOME = REF + '/genome'
DUMP = REF + f'/gff-{MIRBASE_GFF_ID}/unitas-20-11-2026'  # one database per miRNA GFF (#440) and UNITAS dump
DB = DUMP + '/db/genomes/gcf'
CLASSES = ('structural_RNA', 'miRNA', 'hairpin', 'exon', 'TE', 'intron', 'piRNA_cluster')
INDEX = ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt')


def test_route_plans_database_and_runs_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,tesmall', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule tesmall_quant:') == 4
    genome = rule_block(output, 'tesmall_genome')
    assert 'input: data/ext/ensembl/release-116/mus_musculus/GRCm39/fasta/genome.fa' in genome  # the mirflowz reference genome
    assert re.search(rf'bowtie-build -q --threads \d+ {GENOME}/genome.fa {GENOME}/bowtie_index/genome', genome)
    assert f'ln -sfr' in rule_block(output, 'tesmall_genome_links') and f'link={DB}/sequence/${{target#{GENOME}/}}' in output
    db = rule_block(output, 'tesmall_db')
    assert f'--chr-map data/ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus/gff-{MIRBASE_GFF_ID}/chr_map.tsv' in db
    assert '--pirna-clusters data/ext/unitas/20-11-2026/UNITAS_refdump_mus_musculus/piRclDB' in db
    assert f'--trna-archive {REF}/sources/mm39-tRNAs.tar.gz --rmsk {REF}/sources/rmsk.txt.gz' in db
    assert f'--genome-index {GENOME}/bowtie_index/genome' in db and f'--record {DUMP}/tesmall_reference.tsv' in db
    assert 'https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz' in rule_block(output, 'tesmall_gtrnadb')
    assert 'goldenPath/mm39/database/rmsk.txt.gz' in rule_block(output, 'tesmall_rmsk')
    assert 'id=BK000964.3' in rule_block(output, 'tesmall_rdna')
    run = rule_block(output, 'tesmall_quant')
    assert 'TEsmall --no-trim -f "$fastq" -l SR0' in run and f'-g gcf --dbfolder "$db"' in run
    assert '-m 16 -M 36 --mismatch 0 --maxaln 100' in run  # the kit's filtered minimum (NEXTflex v4 fixture)
    assert '-o miRNA hairpin structural_RNA exon TE intron piRNA_cluster' in run  # miRNA first by default
    assert 'input: data/tmp/smallrna/filter/cutadapt/SR0' in run
    assert f'db=$(realpath {DUMP}/db)' in run
    for path in (f'{DB}/annotation/piRNA_cluster.bed', f'{DB}/sequence/bowtie_index/tDNA.rev.2.ebwt', f'{DB}/sequence/genome.fa.fai'):
        assert path in run.split('output:', 1)[0]  # every file TEsmall checks at startup is an input
    returncode, output = dry_run(tmp_path, 'unitas,tesmall', target='bfq_all')
    assert returncode == 0, output
    assert 'tesmall' not in output


def test_each_unitas_dump_has_its_own_database_on_the_shared_genome(tmp_path):
    release = lambda s: s.setdefault('db', {}).setdefault('unitas', {}).update(release='01-01-2027')
    returncode, output = dry_run(tmp_path, 'tesmall', release, printshell=True)
    assert returncode == 0, output
    db = rule_block(output, 'tesmall_db')
    assert '--pirna-clusters data/ext/unitas/01-01-2027/UNITAS_refdump_mus_musculus/piRclDB' in db
    assert f'--record {REF}/gff-{MIRBASE_GFF_ID}/unitas-01-01-2027/tesmall_reference.tsv' in db and f'{DUMP}/' not in output
    assert f'--genome-index {GENOME}/bowtie_index/genome' in db


def test_a_partial_database_is_completed_before_tesmall_runs(tmp_path):
    """With the genome, its links and the database in place except one BED, quantification first rebuilds the database."""
    configure(tmp_path, 'tesmall')
    present = ([f'{GENOME}/genome.fa', f'{GENOME}/genome.fa.fai'] + [f'{GENOME}/bowtie_index/genome.{s}' for s in INDEX] +
               [f'{DB}/sequence/genome.fa', f'{DB}/sequence/genome.fa.fai'] + [f'{DB}/sequence/bowtie_index/genome.{s}' for s in INDEX] +
               [f'{DB}/annotation/{c}.bed' for c in CLASSES if c != 'TE'] +
               [f'{DB}/sequence/{n}.fa{e}' for n in ('rDNA', 'tDNA') for e in ('', '.fai')] +
               [f'{DB}/sequence/bowtie_index/{n}.{s}' for n in ('rDNA', 'tDNA') for s in INDEX] + [f'{DUMP}/tesmall_reference.tsv'])
    for name in present:
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).touch()
    returncode, output = dry_run(tmp_path, 'tesmall')
    assert returncode == 0, output
    assert 'rule tesmall_db:' in output and f'{DB}/annotation/TE.bed' in output


def test_no_shared_minimum_length_gives_tesmall_none(tmp_path):
    fastp = lambda s: s.update(filter={'trim': {'quantifier': 'fastp', 'fastp': {'args': '-L'}}})
    returncode, output = dry_run(tmp_path, 'tesmall', fastp, printshell=True)
    assert returncode == 0, output
    assert '-m 1 -M 36 ' in rule_block(output, 'tesmall_quant')


def test_rat_database_uses_rnor6_sources(tmp_path):
    returncode, output = dry_run(tmp_path, 'tesmall', lambda s: s.update(organism='rattus_norvegicus'), printshell=True)
    assert returncode == 0, output
    assert 'Rnorv6/rn6-tRNAs.tar.gz' in rule_block(output, 'tesmall_gtrnadb') and 'id=NR_046239.1' in rule_block(output, 'tesmall_rdna')
    assert 'data/ext/tesmall/mirbase-23.1_ensembl-104/rattus_norvegicus/' in rule_block(output, 'tesmall_db')


def test_guards_stop_at_load(tmp_path):
    returncode, output = dry_run(tmp_path, 'tesmall', lambda s: s.update(organism='danio_rerio'))
    assert returncode != 0 and 'organism' in output
    returncode, output = dry_run(tmp_path / 'b', 'tesmall', lambda s: s['quant'].setdefault('tesmall', {}).update(max_length=10))
    assert returncode != 0 and 'quant.tesmall.max_length 10 is below the minimum read length' in output
    returncode, output = dry_run(tmp_path / 'c', 'tesmall', lambda s: s['quant'].setdefault('tesmall', {}).update(order='miRNA,hairpin'))
    assert returncode != 0 and 'quant.tesmall.order must list each TEsmall class once' in output
    # TEsmall gets the resolved database path: a symbolic link to a folder named genomes... is caught at load.
    (tmp_path / 'genomes-store').mkdir()
    (tmp_path / 'ext-link').symlink_to(tmp_path / 'genomes-store')
    returncode, output = dry_run(tmp_path / 'd', 'tesmall', lambda s: s.update(ext_dir=str(tmp_path / 'ext-link')))
    assert returncode != 0 and 'contains "genomes", which TEsmall strips from it' in output


def test_tesmall_own_order_can_be_set(tmp_path):
    order = 'structural_RNA,miRNA,hairpin,exon,TE,intron,piRNA_cluster'
    returncode, output = dry_run(tmp_path, 'tesmall', lambda s: s['quant'].setdefault('tesmall', {}).update(order=order), printshell=True)
    assert returncode == 0, output
    assert '-o structural_RNA miRNA hairpin exon TE intron piRNA_cluster' in rule_block(output, 'tesmall_quant')


def test_dotted_sample_names_stop_at_load(tmp_path):
    def rename(settings):
        settings['samples'] = {'S.1' if k == 'SR01' else k: v for k, v in settings['samples'].items()}
    returncode, output = dry_run(tmp_path, 'tesmall', rename)
    assert returncode != 0 and 'TEsmall cuts sample labels at dots' in output


def test_counts_join_samples_and_reject_unexpected_tables(tmp_path):
    (tmp_path / 'a.txt').write_text('fid\tftype\tS1\nmmu-miR-1a-3p\tmiRNA\t5\nL1:L1:LINE\tsense_TE\t2\n')
    (tmp_path / 'b.txt').write_text('fid\tftype\tS2\nmmu-miR-1a-3p\tmiRNA\t7\nmmu-mir-1a-1\thairpin\t1\n')
    script = [sys.executable, str(ROOT / 'smallrna/rules/quant/scripts/tesmall_counts.py')]
    result = subprocess.run(script + ['--samples', 'S1', 'S2', '--tables', str(tmp_path / 'a.txt'), str(tmp_path / 'b.txt')],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ['ftype\tfid\tS1\tS2', 'hairpin\tmmu-mir-1a-1\t0\t1', 'miRNA\tmmu-miR-1a-3p\t5\t7',
                                          'sense_TE\tL1:L1:LINE\t2\t0']
    result = subprocess.run(script + ['--samples', 'S2', '--tables', str(tmp_path / 'a.txt')], capture_output=True, text=True)
    assert result.returncode != 0 and 'expected columns fid, ftype, S2' in result.stderr


def test_tesmall_runs_with_its_temporary_files_in_the_work_directory_and_is_checked(tmp_path):
    """#419: TMPDIR inside the sample's work directory, and tesmall_check.py after TEsmall."""
    returncode, output = dry_run(tmp_path, 'tesmall', printshell=True)
    assert returncode == 0, output
    run = rule_block(output, 'tesmall_quant')
    assert 'cd data/tmp/smallrna/quant/tesmall/SR0' in run and ' && rm -rf tmp && mkdir tmp && TMPDIR=$PWD/tmp ' in run
    assert '--verbose 2 > "$log" 2>&1 && rm -rf tmp && python "$script" --work . --sample SR0' in run
    assert re.search(r'script=\$\(realpath \S+/smallrna/rules/quant/scripts/tesmall_check.py\)', run)
    assert 'fastq=$(realpath data/tmp/smallrna/filter/cutadapt/SR0' in run  # :q quoting leaves plain paths as they are


CHECK = ROOT / 'smallrna/rules/quant/scripts/tesmall_check.py'


def bowtie_log(aligned=60):
    return (f'# reads processed: 100\n# reads with at least one alignment: {aligned} ({aligned}.00%)\n'
            f'# reads that failed to align: {100 - aligned} ({100 - aligned}.00%)\nReported {aligned} alignments\n')


def tesmall_work(work, sample='S1', anno_rows=3, aligned_rows=5, genome_aligned=60, log_lines=('INFO TEsmall run completed',)):
    """A TEsmall work directory as a finished run leaves it, reduced to what tesmall_check.py reads."""
    work.mkdir(parents=True, exist_ok=True)
    (work / 'tesmall.log').write_text('\n'.join(log_lines) + '\n')
    for step in ('rRNA', 'genome', 'tRNA'):
        (work / f'{sample}.{step}.log').write_text(bowtie_log(genome_aligned if step == 'genome' else 60))
    (work / f'{sample}.aligned.rinfo').write_text('id\tlength\tcount\tmapper\n' + 'r\t22\t1\tunique\n' * aligned_rows)
    (work / f'{sample}.anno').write_text('rid\trchr\n' + 'r\t1\n' * anno_rows)
    return work


def check(work, sample='S1'):
    return subprocess.run([sys.executable, str(CHECK), '--work', str(work), '--sample', sample, '--log', str(work / 'tesmall.log')],
                          capture_output=True, text=True)


def test_check_passes_a_complete_run_and_one_without_aligned_reads(tmp_path):
    assert check(tesmall_work(tmp_path / 'a')).returncode == 0  # Bowtie's "failed to align" summary is not an error
    assert check(tesmall_work(tmp_path / 'b', anno_rows=0, aligned_rows=0, genome_aligned=0)).returncode == 0
    work = tesmall_work(tmp_path / 'c')  # a nonfatal Bowtie warning mentioning "error" in a file name
    (work / 'S1.genome.log').write_text('Warning: skipping read in error_reads.fq\n' + bowtie_log())
    assert check(work).returncode == 0


@pytest.mark.parametrize('change, message', [
    (lambda w: tesmall_work(w, anno_rows=0), '60 reads aligned to the genome but S1.anno is empty'),
    (lambda w: (tesmall_work(w), (w / 'S1.aligned.rinfo').unlink()), '60 reads aligned to the genome but S1.aligned.rinfo is missing'),
    (lambda w: tesmall_work(w, log_lines=('sort: read failed: /x/tmp/sortQ: Input/output error', 'INFO TEsmall run completed')),
     'tesmall.log: sort: read failed'),  # a partial, non-empty annotation passes the row check; the log does not
    (lambda w: (tesmall_work(w), (w / 'S1.genome.log').write_text('# reads processed: 100\n')), 'S1.genome.log: no aligned-read count'),
    (lambda w: tesmall_work(w, log_lines=('INFO Annotating...', 'sort: write failed: /tmp/sortX: Disk quota exceeded',
                                          'INFO TEsmall run completed')), 'tesmall.log: sort: write failed'),
    (lambda w: tesmall_work(w, log_lines=('[E::hts_open_format] Failed to open file "S1.genome.bam"',)), '[E::hts_open_format]'),
    (lambda w: (tesmall_work(w), (w / 'S1.rRNA.log').write_text('Could not locate a Bowtie index corresponding to basename "x"\n')),
     'S1.rRNA.log: Could not locate a Bowtie index'),
    (lambda w: (tesmall_work(w), (w / 'S1.genome.log').write_text('')), 'S1.genome.log: no "# reads processed" summary'),
    (lambda w: (tesmall_work(w), (w / 'S1.tRNA.log').unlink()), 'S1.tRNA.log: missing'),
    (lambda w: (tesmall_work(w), (w / 'tesmall.log').write_text('')), 'missing or empty'),
])
def test_check_stops_a_run_whose_steps_failed(tmp_path, change, message):
    change(tmp_path / 'w')
    result = check(tmp_path / 'w')
    assert result.returncode != 0 and message in result.stderr, result.stderr


def test_a_newer_mirna_gff_or_gtf_does_not_rebuild_the_shared_genome_indexes(tmp_path):
    """tesmall_genome reads only the genome (#433), so a new miRNA GFF3 or GTF rebuilds what is built from them, not
    the genome FASTA, its Bowtie index or the biocore_htseq Bowtie 2 index built from that FASTA."""
    configure(tmp_path, 'tesmall,biocore_htseq')
    targets = ('tesmall_genome', 'biocore_htseq_bowtie2_index', 'tinyrna_gff', 'tesmall_db')
    returncode, output = snakemake(tmp_path, '-n', *targets)
    assert returncode == 0, output
    built = {path for line in re.findall(r'^    (?:output|log|benchmark): (.+)$', output, re.M) for path in line.split(', ')}
    gtf = 'data/ext/ensembl/release-116/mus_musculus/GRCm39/anno/mus_musculus.GRCm39.116.gtf'
    assert {MIRBASE_GFF, gtf} <= built
    earlier = time.time() - 86400
    for name in sorted(built):  # every file the targets need (a missing benchmark file reruns its job), all built at once
        path = tmp_path / name
        if any(other.startswith(name + '/') for other in built):
            path.mkdir(parents=True, exist_ok=True)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
    for name in built:
        os.utime(tmp_path / name, (earlier, earlier))
    for name in (MIRBASE_GFF, gtf):  # a new annotation
        os.utime(tmp_path / name, (earlier + 3600, earlier + 3600))
    returncode, output = snakemake(tmp_path, '-n', *targets)
    assert returncode == 0, output
    planned = set(re.findall(r'^(?:local)?rule (\w+):', output, re.M))
    assert {'tesmall_db', 'mirflowz_chr_map', 'tinyrna_gff'} <= planned, output
    assert not planned & {'tesmall_genome', 'biocore_htseq_bowtie2_index'}, output


def planned_outputs(output):
    """{rule: {output paths}} of a dry run's jobs."""
    planned = {}
    for block in re.split(r'^(?:local)?rule ', output, flags=re.M)[1:]:
        match = re.search(r'^    output: (.+)$', block, re.M)
        if match:
            planned.setdefault(block.split(':', 1)[0], set()).update(match.group(1).split(', '))
    return planned


def test_two_correction_sets_in_one_ext_dir_each_use_their_own_derived_files(tmp_path):
    """#440: two checkouts whose miRBase GFF3 corrections differ share one EXT_DIR. Everything built from the corrected
    GFF3 is named by its identity, so each plans and keeps its own; what is built from the genome alone is shared."""
    other = tmp_path / 'other'  # a second checkout, with one more byte in the correction table
    for name in subprocess.run(['git', 'ls-files'], cwd=ROOT, capture_output=True, text=True, check=True).stdout.splitlines():
        (other / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, other / name)
    table = other / 'gcfdb/tables/mirbase_gff_corrections.tsv'
    table.write_bytes(table.read_bytes() + b'\n')
    ids = {'a': MIRBASE_GFF_ID, 'b': mirbase_gff_id(other)}
    assert ids['a'] != ids['b']
    ext = tmp_path / 'ext'
    derived = ('mirbase_genes', 'mirflowz_chr_map', 'mirflowz_annotation', 'tinyrna_gff', 'tesmall_db', 'biocore_miraligner_db')
    genome = ('mirflowz_prepare', 'tesmall_genome', 'biocore_htseq_bowtie2_index')
    for name, root in (('a', ROOT), ('b', other)):
        configure(tmp_path / name, 'mirflowz,tesmall,tinyrna,biocore_htseq,biocore_miraligner', lambda s: s.update(ext_dir=str(ext)))
        (tmp_path / name / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                                   f'include: {str(root / "smallrna/smallrna.smk")!r}\n')
    plans = {}
    for name in ids:
        returncode, output = snakemake(tmp_path / name, '-n', *derived, *genome)
        assert returncode == 0, output
        plans[name] = planned_outputs(output)
    for rule in derived:  # each checkout's own files, named by its corrections
        for name, gff_id in ids.items():
            assert all(f'gff-{gff_id}/' in path or f'genes.{gff_id}.' in path for path in plans[name][rule]), (rule, plans[name][rule])
    for rule in genome:  # one copy for both
        assert plans['a'][rule] == plans['b'][rule] and not any('gff-' in path for path in plans['a'][rule]), rule

    def build(name, when):
        returncode, output = snakemake(tmp_path / name, '-n', *derived, *genome)
        assert returncode == 0, output
        made = {p for line in re.findall(r'^    (?:output|log|benchmark): (.+)$', output, re.M) for p in line.split(', ')}
        for path in map(Path, made):
            if any(other.startswith(str(path) + '/') for other in made) or path.name in ('prepared', 'annotation'):
                path.mkdir(parents=True, exist_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
            os.utime(path, (when, when))
        return set(re.findall(r'^(?:local)?rule (\w+):', output, re.M))

    earlier = time.time() - 86400
    assert set(derived) | set(genome) <= build('a', earlier)
    returncode, output = snakemake(tmp_path / 'a', '-n', *derived, *genome)
    assert returncode == 0 and 'Nothing to be done' in output, output
    planned = build('b', earlier + 3600)  # the other corrections, built later in the same EXT_DIR
    assert set(derived) <= planned and not planned & set(genome), planned
    returncode, output = snakemake(tmp_path / 'a', '-n', *derived, *genome)
    assert returncode == 0 and 'Nothing to be done' in output, output  # a keeps its own files, older than b's


RNACENTRAL = DUMP + '-rnacentral-28.0'  # #378: the database with RNAcentral's structural entries, beside the base one


def test_rnacentral_structural_entries_make_a_database_of_their_own(tmp_path):
    """quant.tesmall.structural_rnacentral (#378): RNAcentral's structural entries in a database derived from the base
    one, which TEsmall then reads; the base database is built as before, and without the setting nothing of it is
    planned. Rat (Rnor_6.0) stops at load."""
    on = lambda s: s['quant'].setdefault('tesmall', {}).update(structural_rnacentral=True)
    returncode, output = dry_run(tmp_path, 'tesmall', on, printshell=True)
    assert returncode == 0, output
    assert ('https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/mus_musculus.GRCm39.gff3.gz'
            in rule_block(output, 'tesmall_rnacentral'))
    derived = rule_block(output, 'tesmall_rnacentral_db')
    assert f'--base {DB} ' in derived and f'--out {RNACENTRAL}/db/genomes/gcf ' in derived
    assert f'--record {RNACENTRAL}/tesmall_reference.tsv' in derived and f'--base-record {DUMP}/tesmall_reference.tsv' in derived
    assert f'--rnacentral {RNACENTRAL}/sources/mus_musculus.GRCm39.gff3.gz' in derived   # in its own folder, not the shared sources/
    for path in ('annotation/structural_RNA.bed', 'annotation/piRNA_cluster.bed', 'sequence/genome.fa', 'sequence/bowtie_index/tDNA.rev.2.ebwt'):
        assert f'{RNACENTRAL}/db/genomes/gcf/{path}' in derived.split('output:', 1)[1].split('shell:', 1)[0]
    run = rule_block(output, 'tesmall_quant')
    assert f'db=$(realpath {RNACENTRAL}/db)' in run and f'{RNACENTRAL}/db/genomes/gcf/annotation/structural_RNA.bed' in run.split('output:', 1)[0]
    assert 'rule tesmall_db:' in output   # the base database it is derived from
    returncode, output = dry_run(tmp_path / 'off', 'tesmall', printshell=True)
    assert returncode == 0, output
    assert 'rule tesmall_rnacentral' not in output and f'db=$(realpath {DUMP}/db)' in rule_block(output, 'tesmall_quant')
    returncode, output = dry_run(tmp_path / 'rat', 'tesmall', lambda s: (on(s), s.update(organism='rattus_norvegicus')))
    assert returncode != 0 and 'RNAcentral 28.0 has coordinates for homo_sapiens and mus_musculus on GRCh38 and GRCm39' in output
    returncode, output = dry_run(tmp_path / 'text', 'tesmall', lambda s: s['quant'].setdefault('tesmall', {}).update(structural_rnacentral='no'))
    assert returncode != 0 and "quant.tesmall.structural_rnacentral must be true or false, got 'no'" in output
