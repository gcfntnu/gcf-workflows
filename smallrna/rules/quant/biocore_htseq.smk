#-*- mode:snakemake -*-
"""The genome counter of an earlier in-house small-RNA pipeline, Bowtie 2 and htseq-count, as an optional quant_all
route (quant.method: biocore_htseq, alone or combined, e.g. unitas,biocore_htseq; #425). Delivery (bfq_all, multiqc_report) does not
use it.

The pipeline's own trimming is replaced by the shared filtered reads, as for the other routes. Per sample:

1. the reads are collapsed to unique sequences with fastx_collapser (names `<rank>-<reads>`);
2. sequences of at least 15 nt are aligned end to end to the genome with Bowtie 2, `-k 10 -f` with default
   sensitivity: up to 10 alignments per sequence, mismatches within Bowtie 2's score threshold. With
   quant.biocore_htseq.alignment `local` (#417), Bowtie 2 aligns locally (`--local`, minimum score 2 per base for
   `local_min_aligned` of the read, the match bonus, so that share must align), and alignments whose read 5' end is
   soft-clipped by more than `local_max_5p_clip` bases or that align less than `local_min_aligned` of the read are
   dropped before annotation: a clipped 3' tail passes step 4's perfect-alignment test, which reads the aligned part
   only. Of the alignments left, a read keeps only those with the highest score (AS), so an exact read is not also
   counted at a relative's locus that it matches clipped. The default local_min_aligned of 0.9 covers tails of up to
   two bases on reads of 20-29 nt;
3. htseq-count annotates every alignment (`-a 0 -s yes -t miRNA -i Name`: union overlap mode, stranded, MAPQ 0
   accepted), writing the feature into the record's XF tag. Bowtie 2 writes no NH tag, so each reported alignment is
   annotated on its own. Secondary (and supplementary) alignments are scored, as HTSeq 0.9.1, the pipeline's version,
   did by default; HTSeq 2 ignores them by default, which would leave out every alignment of a read but its primary
   one, so the route sets both explicitly;
4. htseq_counts.py keeps the alignments the pipeline took as perfect (no mismatch in the MD tag's first run, see the
   script) and counts each (read, miRNA name) pair once, ambiguous feature sets expanded: a read counts in full for
   every miRNA name it overlaps at any of its loci, never divided.

Per sample `mature.tsv` (name, count), `isomir.tsv` (name_sequence, count) and `stats.tsv`; `htseq_counts.tsv` joins
the samples' mature counts. The miRNA annotation is tiny-count's renamed GFF (gcfdb/tinyrna.db) without its rows that
have no strand, which htseq-count's stranded mode cannot place and stops on; the corrected miRBase GFF3 has none
(gcfdb/mirbase.db, #429), a custom reference (`db.mirflowz.custom`) may.

ncRNA tables (#503). The pipeline's annotation added RNAcentral transcripts and Ensembl genes to the miRNAs, all as
feature type `miRNA`, and one htseq-count pass annotated the alignments with all three; where gcfdb/biocore_htseq.db
builds that annotation (`BIOCORE_HTSEQ_CLASSES`: human on the Ensembl reference, with the pipeline's RNAcentral 9.0),
step 3 uses it. In union mode the extra features only turn a miRNA hit into `__ambiguous[...]`, which step 4 expands
back, so the miRNA tables do not depend on them (#427 and #503 showed identical tables). From the same records htseq_counts.py `ncrna` writes the pipeline's
ncRNA tables per sample: `ncrna_types.tsv` (RNAcentral types and Ensembl biotypes, each read divided evenly between its
distinct types), `ncrna.tsv` (URS/ENS feature and sequence) and `ncrna_stats.tsv` (RNACentral, tRNAs, rRNAs, snoRNAs,
snRNAs, lncRNAs, by substring as the pipeline's stats); `ncrna_types.tsv` in the route's folder joins the samples'
classes. The pipeline's quirks are kept: the split between the miRNA and ncRNA paths is per alignment record, so a
read with alignments on both counts in both; RNAcentral types and Ensembl biotypes share one class list. Elsewhere the
route annotates with the miRNA GFF alone (`mirna.gff3`) and writes no ncRNA tables. The genome's Bowtie 2 index is
gcfdb/biocore_htseq.db.
"""

include:
    join(GCFDB_DIR, 'biocore_htseq.db')

BIOCORE_HTSEQ_OUT = join(QUANT_INTERIM, 'biocore_htseq')
BIOCORE_HTSEQ_MIN_LENGTH = 15  # collapsed sequences shorter than this are not aligned
BIOCORE_HTSEQ_BOWTIE2 = '-k 10 -f'
BIOCORE_HTSEQ_COUNT = '-f sam -a 0 -s yes -t miRNA -i Name --secondary-alignments score --supplementary-alignments score'
BIOCORE_HTSEQ_CONF = config['quant'].get('biocore_htseq') or {}
if not isinstance(BIOCORE_HTSEQ_CONF, dict):
    raise ValueError('quant.biocore_htseq must be a mapping of settings, got {!r}'.format(BIOCORE_HTSEQ_CONF))


def biocore_htseq_local():
    """(Bowtie 2 options, filter) of quant.biocore_htseq.alignment (#417), joined into the rule's command here so that
    the default (end_to_end) keeps it as before: for `local`, --local with a minimum score of 2 (the match bonus) per
    base of local_min_aligned of the read, and an awk filter of the SAM records, unmapped ones kept, that drops
    alignments clipping more than local_max_5p_clip bases at the read's 5' end (the CIGAR's start on the forward strand,
    its end on the reverse) or aligning less than local_min_aligned of the read, and of the rest keeps a read's
    alignments with the highest score (AS; a read's records are consecutive)."""
    mode = BIOCORE_HTSEQ_CONF.get('alignment', 'end_to_end')
    if mode == 'end_to_end':
        return '', ''
    if mode != 'local':
        raise ValueError('quant.biocore_htseq.alignment must be end_to_end or local, got {!r}'.format(mode))
    share, clip = BIOCORE_HTSEQ_CONF.get('local_min_aligned'), BIOCORE_HTSEQ_CONF.get('local_max_5p_clip')
    if isinstance(share, bool) or not isinstance(share, (int, float)) or not 0.5 <= share <= 1:
        raise ValueError('quant.biocore_htseq.local_min_aligned must be a number from 0.5 to 1, got {!r}'.format(share))
    if isinstance(clip, bool) or not isinstance(clip, int) or not 0 <= clip <= 5:
        raise ValueError('quant.biocore_htseq.local_max_5p_clip must be a whole number from 0 to 5, got {!r}'.format(clip))
    awk = ("awk -F '\\t' -v OFS='\\t' -v share={share} -v clip={clip} "
           "'function flush(  j) {{{{ for (j = 1; j <= k; j++) print rec[j]; k = 0 }}}} "
           "/^@/ {{{{ print; next }}}} "
           "int($2 / 4) % 2 {{{{ flush(); name = \"\"; print; next }}}} "
           "{{{{ if ($1 != name) {{{{ flush(); name = $1; best = -1e9 }}}} "
           "n = length($10); lead = ($6 ~ /^[0-9]+S/) ? $6 + 0 : 0; tail = 0; "
           "if (match($6, /[0-9]+S$/)) tail = substr($6, RSTART, RLENGTH - 1) + 0; "
           "five = int($2 / 16) % 2 ? tail : lead; "
           "if (five > clip || n - lead - tail < share * n - 1e-9) next; "
           "as = -1e9; for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) as = substr($i, 6) + 0; "
           "if (as > best) {{{{ best = as; k = 0 }}}} "
           "if (as == best) rec[++k] = $0 }}}} "
           "END {{{{ flush() }}}}'").format(share=share, clip=clip)
    # a seed of 16 (default 20 in local mode), so a 20-nt read with a 2-base tail still has an exact seed
    return ' --local -L 16 --score-min L,0,{:g}'.format(2 * share), ' | ' + awk


BIOCORE_HTSEQ_LOCAL, BIOCORE_HTSEQ_FILTER = biocore_htseq_local()  # also when its outputs are requested without quant.method selecting it
BIOCORE_HTSEQ_LOCAL_MODE = BIOCORE_HTSEQ_CONF.get('alignment', 'end_to_end') == 'local'


def biocore_htseq_settings():
    """The settings that decide the biocore_htseq route's counting, recorded with its table (quant_settings, #417)."""
    return [('minimum_length', str(BIOCORE_HTSEQ_MIN_LENGTH), 'fixed'),
            ('alignment', str(BIOCORE_HTSEQ_CONF.get('alignment', 'end_to_end')), 'quant.biocore_htseq.alignment'),
            # recorded as unused in end_to_end mode, so no benchmark variant can vary them alone
            ('local_min_aligned', str(BIOCORE_HTSEQ_CONF.get('local_min_aligned')) if BIOCORE_HTSEQ_LOCAL_MODE else 'unused (end_to_end)',
             'quant.biocore_htseq.local_min_aligned'),
            ('local_max_5p_clip', str(BIOCORE_HTSEQ_CONF.get('local_max_5p_clip')) if BIOCORE_HTSEQ_LOCAL_MODE else 'unused (end_to_end)',
             'quant.biocore_htseq.local_max_5p_clip'),
            ('bowtie2', BIOCORE_HTSEQ_BOWTIE2 + BIOCORE_HTSEQ_LOCAL, 'quant.biocore_htseq.alignment and fixed'),
            ('htseq_count', BIOCORE_HTSEQ_COUNT, 'fixed'),
            ('counting', 'perfect alignments only (MD tag); a read counts in full for every miRNA it overlaps (htseq_counts.py)', 'fixed'),
            ('ncrna_annotation', 'RNAcentral {} transcripts and the reference\'s Ensembl genes (#503)'.format(BIOCORE_HTSEQ_RNACENTRAL_RELEASE)
             if BIOCORE_HTSEQ_CLASSES else 'none (miRNA GFF only)', 'fixed per organism and reference')]


rule biocore_htseq_gff:
    input:
        rules.tinyrna_gff.output.gff
    output:
        join(BIOCORE_HTSEQ_OUT, 'mirna.gff3')
    params:
        script = src_gcf('scripts/htseq_counts.py')
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_gff.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} gff --input {input:q} > {output:q}'


rule biocore_htseq_collapse:
    input:
        fastq = lambda wildcards: get_filtered_fastq(wildcards)['R1']
    output:
        join(BIOCORE_HTSEQ_OUT, '{sample}', 'collapsed.fa')
    log:
        join(BIOCORE_HTSEQ_OUT, '{sample}', 'collapse.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_collapse', '{sample}.tsv')
    container:
        'docker://' + config['docker']['fastx_toolkit']
    shell:
        # fastx_collapser stops on an empty FASTQ; a sample without filtered reads has no sequences (and zero counts).
        'if [ -s {input.fastq:q} ]; then fastx_collapser -Q 33 -i {input.fastq:q} -o {output:q} 2> {log:q}; '
        'else : > {output:q}; echo "no filtered reads" > {log:q}; fi'


rule biocore_htseq_align:
    input:
        reads = rules.biocore_htseq_collapse.output,
        index = rules.biocore_htseq_bowtie2_index.output
    output:
        join(BIOCORE_HTSEQ_OUT, '{sample}', 'alignments.sam.gz')
    params:
        index = BIOCORE_HTSEQ_INDEX
    threads:
        4
    resources:
        # Bowtie 2 holds the genome index in memory, about 3.5 GB for the human genome.
        mem_mb = 5000
    log:
        join(BIOCORE_HTSEQ_OUT, '{sample}', 'bowtie2.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_align', '{sample}.tsv')
    container:
        'docker://' + config['docker']['biocore_htseq_bowtie2']
    shell:
        # fastx_collapser writes each sequence on one line after its name.
        "awk 'NR % 2 == 1 {{name = $0; next}} length($0) >= " + str(BIOCORE_HTSEQ_MIN_LENGTH) + " {{print name; print}}' {input.reads:q} | "
        'bowtie2 -p {threads} ' + BIOCORE_HTSEQ_BOWTIE2 + BIOCORE_HTSEQ_LOCAL + ' -x {params.index:q} -U - 2> {log:q}' + BIOCORE_HTSEQ_FILTER +
        ' | gzip > {output:q}'


rule biocore_htseq_annotate:
    input:
        alignments = rules.biocore_htseq_align.output,
        # the pipeline's annotation with the RNAcentral and Ensembl layers (#503) where it is built, else the miRNA GFF
        gff = BIOCORE_HTSEQ_ANNOTATION if BIOCORE_HTSEQ_CLASSES else rules.biocore_htseq_gff.output
    output:
        join(BIOCORE_HTSEQ_OUT, '{sample}', 'annotated.sam.gz')
    params:
        samout = join(BIOCORE_HTSEQ_OUT, '{sample}', 'annotated.sam')
    resources:
        mem_mb = 1000  # about 140 MB with the ncRNA annotation (#503)
    log:
        join(BIOCORE_HTSEQ_OUT, '{sample}', 'htseq-count.log')
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_annotate', '{sample}.tsv')
    container:
        'docker://' + config['docker']['htseq']
    shell:
        # The uncompressed samout is not a declared output: removed if htseq-count fails, compressed if it succeeds.
        '(gzip -dc {input.alignments:q} | '
        'htseq-count ' + BIOCORE_HTSEQ_COUNT + ' '
        '--samout {params.samout:q} - {input.gff:q} > {log:q} 2>&1) || {{ rm -f {params.samout:q}; exit 1; }}; '
        'gzip -f {params.samout:q}'


rule biocore_htseq_count:
    input:
        collapsed = rules.biocore_htseq_collapse.output,
        annotated = rules.biocore_htseq_annotate.output
    output:
        mature = join(BIOCORE_HTSEQ_OUT, '{sample}', 'mature.tsv'),
        isomir = join(BIOCORE_HTSEQ_OUT, '{sample}', 'isomir.tsv'),
        stats = join(BIOCORE_HTSEQ_OUT, '{sample}', 'stats.tsv')
    params:
        script = src_gcf('scripts/htseq_counts.py'),
        species = MIR_ORG.get(ORG, '')
    resources:
        mem_mb = 1000
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_count', '{sample}.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} count --species {params.species:q} --collapsed {input.collapsed:q} '
        '--annotated {input.annotated:q} --mature {output.mature:q} --isomir {output.isomir:q} --stats {output.stats:q}'


rule biocore_htseq_counts:
    input:
        expand(rules.biocore_htseq_count.output.mature, sample=SAMPLES)
    output:
        join(BIOCORE_HTSEQ_OUT, 'htseq_counts.tsv')
    params:
        script = src_gcf('scripts/htseq_counts.py'),
        samples = SAMPLES
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_counts.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} join --samples {params.samples:q} --tables {input:q} > {output:q}'


rule biocore_htseq_ncrna:
    input:
        rules.biocore_htseq_annotate.output
    output:
        types = join(BIOCORE_HTSEQ_OUT, '{sample}', 'ncrna_types.tsv'),
        features = join(BIOCORE_HTSEQ_OUT, '{sample}', 'ncrna.tsv'),
        stats = join(BIOCORE_HTSEQ_OUT, '{sample}', 'ncrna_stats.tsv')
    params:
        script = src_gcf('scripts/htseq_counts.py'),
        species = MIR_ORG.get(ORG, '')
    resources:
        mem_mb = 1000  # about 50 MB on a plasma sample
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_ncrna', '{sample}.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} ncrna --species {params.species:q} --annotated {input:q} --types {output.types:q} '
        '--features {output.features:q} --stats {output.stats:q}'


rule biocore_htseq_ncrna_types:
    input:
        expand(rules.biocore_htseq_ncrna.output.types, sample=SAMPLES)
    output:
        join(BIOCORE_HTSEQ_OUT, 'ncrna_types.tsv')
    params:
        script = src_gcf('scripts/htseq_counts.py'),
        samples = SAMPLES
    benchmark:
        join(QUANT_BENCHMARKS, 'biocore_htseq_ncrna_types.tsv')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script:q} join --key type --samples {params.samples:q} --tables {input:q} > {output:q}'


rule biocore_htseq_all:
    input:
        rules.biocore_htseq_counts.output,
        expand(rules.biocore_htseq_count.output, sample=SAMPLES),
        # the ncRNA tables (#503) where the annotation has the RNAcentral and Ensembl layers
        rules.biocore_htseq_ncrna_types.output if BIOCORE_HTSEQ_CLASSES else [],
        expand(rules.biocore_htseq_ncrna.output, sample=SAMPLES) if BIOCORE_HTSEQ_CLASSES else []
