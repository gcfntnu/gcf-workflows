#-*- mode: snakemake -*-
"""miRTrace contamination check and miRNA complexity on the filtered reads (#352).

miRTrace reads the reads the quantifiers read (filtered_fastq_pattern: the fastp or cutadapt route, or with
filter.trim.quantifier: skip the merged reads unchanged) without an adapter and with protocol illumina, so it trims
nothing: the filter has already removed adapters and random ends. Its fixed filters (reads under 18 nt, of low
complexity or low PHRED) cannot be switched off; they only leave reads out of miRTrace's own statistics, not out of
the counts. The MultiQC report shows its contamination check and miRNA complexity only
(misc/multiqc/multiqc_config-smallrna.yaml).
"""

MIRTRACE_DIR = join(QC_INTERIM, 'mirtrace')
MIRTRACE_SPECIES = {'homo_sapiens': 'hsa', 'mus_musculus': 'mmu', 'rattus_norvegicus': 'rno'}
MIRTRACE_FASTQ = expand(filtered_fastq_pattern(), sample=SAMPLES)

def mirtrace_title():
    """One title string; several projects are joined with commas, as in the MultiQC report title."""
    project_id = config.get('project_id', 'GCF-0000-00')
    if isinstance(project_id, (list, tuple)):
        return ','.join(map(str, project_id))
    return str(project_id)

def mirtrace_species():
    species = MIRTRACE_SPECIES.get(config['organism'])
    if species is None:
        raise ValueError('miRTrace QC supports {}, got organism {!r}'.format(', '.join(MIRTRACE_SPECIES), config['organism']))
    return species

rule mirtrace_config:
    """One row per sample: filtered FASTQ and name, no adapter column (miRTrace's default: no adapter trimming)."""
    output:
        join(MIRTRACE_DIR, 'mirtrace.config')
    params:
        rows = [[os.path.abspath(fastq), sample] for sample, fastq in zip(SAMPLES, MIRTRACE_FASTQ)]
    run:
        with open(output[0], 'w') as fh:
            for row in params.rows:
                fh.write(','.join(row) + '\n')

rule qc_mirtrace:
    input:
        config = rules.mirtrace_config.output,
        fastq = MIRTRACE_FASTQ
    output:
        results = join(MIRTRACE_DIR, 'mirtrace-results.json'),
        report = join(MIRTRACE_DIR, 'mirtrace-report.html'),
        length = join(MIRTRACE_DIR, 'mirtrace-stats-length.tsv'),
        contamination = join(MIRTRACE_DIR, 'mirtrace-stats-contamination_basic.tsv'),
        complexity = join(MIRTRACE_DIR, 'mirtrace-stats-mirna-complexity.tsv')
    params:
        outdir = MIRTRACE_DIR,
        species = mirtrace_species(),
        title = mirtrace_title()
    threads:
        8
    container:
        'docker://' + config['docker']['mirtrace']
    shell:
        'mirtrace qc --force '
        '-s {params.species} -p illumina -c {input.config} '
        '-o {params.outdir} --title {params.title:q} -t {threads}'
