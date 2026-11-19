#-*- mode:snakemake -*-
"""nf-core/smrnaseq 2.4.1's mirtop counting as an optional quant_all route (quant.method: smrnaseq, alone or combined,
e.g. unitas,smrnaseq; #403). Delivery (bfq_all, multiqc_report) does not use it.

The pipeline's own read preparation (fastp, UMI handling, contaminant filtering) is replaced by the shared filtered
reads, as for the other routes; the development run in #403 showed that with `--skip_fastp` the pipeline counts
exactly these reads. Its MIRNA_QUANT path to `mirna.tsv` is reproduced with the same tools, images and arguments:

- the organism's miRBase hairpins (db.mirbase.release) prepared as the pipeline prepares them (smrnaseq_hairpin.py)
  and indexed with Bowtie 1;
- per sample: reads collapsed to unique sequences (`seqcluster collapse -m 1 --min_size 15`), aligned to the hairpins
  (`bowtie -t -k 50 --best --strata -e 99999 --chunkmbs 2048`: Bowtie's default -n 2 mode, so up to 2 mismatches in
  the first 28 bases, its seed, and with -e 99999 no limit beyond it; the alignments of the best stratum, at most 50),
  annotated (`mirtop gff`) and counted per isomiR (`mirtop counts --add-extra`);
- the samples joined and summed per miRNA (smrnaseq_counts.py): a collapsed read annotated to several miRNAs counts
  in full for the one whose line sorts first, as in the pipeline.

The pipeline's other count table, Bowtie to mature and then hairpin sequences counted with `samtools idxstats`, is not
reproduced: idxstats counts alignment records, so reads matching several mature sequences count more than once, and
isomiRs longer than the mature sequence are missed (#403).
"""

SMRNASEQ_SPECIES = {'homo_sapiens': 'hsa', 'mus_musculus': 'mmu', 'rattus_norvegicus': 'rno'}
SMRNASEQ_SELECTED = 'smrnaseq' in QUANT_SELECTED
SMRNASEQ_OUT = join(QUANT_INTERIM, 'smrnaseq')
SMRNASEQ_REF = join(SMRNASEQ_OUT, 'reference')
SMRNASEQ_INDEX = ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt')

if SMRNASEQ_SELECTED and config['organism'] not in SMRNASEQ_SPECIES:
    raise ValueError('quant.method smrnaseq supports organism {}, got {!r}'.format(', '.join(SMRNASEQ_SPECIES), config['organism']))
SMRNASEQ_SPS = SMRNASEQ_SPECIES.get(config['organism'], '')


rule smrnaseq_hairpin:
    input:
        rules.mirbase_hairpin_fasta.output.fasta
    output:
        join(SMRNASEQ_REF, 'hairpin.fa')
    params:
        script = src_gcf('scripts/smrnaseq_hairpin.py'),
        species = SMRNASEQ_SPS
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} --species {params.species} {input} > {output}'


rule smrnaseq_hairpin_index:
    input:
        rules.smrnaseq_hairpin.output
    output:
        expand(join(SMRNASEQ_REF, 'bowtie', 'hairpin.{suffix}'), suffix=SMRNASEQ_INDEX)
    params:
        prefix = join(SMRNASEQ_REF, 'bowtie', 'hairpin')
    threads:
        4
    log:
        join(SMRNASEQ_REF, 'bowtie', 'hairpin.log')
    container:
        'docker://' + config['docker']['smrnaseq_bowtie']
    shell:
        'bowtie-build --threads {threads} {input} {params.prefix} > {log} 2>&1'


rule smrnaseq_collapse:
    input:
        lambda wildcards: get_filtered_fastq(wildcards)['R1']
    output:
        join(SMRNASEQ_OUT, 'collapsed', '{sample}.fastq')
    params:
        work = join(SMRNASEQ_OUT, 'collapsed', '{sample}.work')
    log:
        join(SMRNASEQ_OUT, 'collapsed', '{sample}.log')
    container:
        'docker://' + config['docker']['smrnaseq_seqcluster']
    shell:
        # seqcluster names its output after the input file (<name>_trimmed.fastq), so each sample gets its own directory.
        'rm -rf "{params.work}" && seqcluster collapse -m 1 --min_size 15 -f "{input}" -o "{params.work}" > "{log}" 2>&1 && '
        'mv "{params.work}"/*_trimmed.fastq "{output}" && rm -r "{params.work}"'


rule smrnaseq_align:
    input:
        fastq = rules.smrnaseq_collapse.output,
        index = rules.smrnaseq_hairpin_index.output
    output:
        # the BAM's name is the sample name in mirtop's tables
        join(SMRNASEQ_OUT, 'bam', '{sample}.bam')
    params:
        index = rules.smrnaseq_hairpin_index.params.prefix,
        args = '-t -k 50 --best --strata -e 99999 --chunkmbs 2048'
    threads:
        4
    log:
        join(SMRNASEQ_OUT, 'bam', '{sample}.log')
    container:
        'docker://' + config['docker']['smrnaseq_bowtie']
    shell:
        # samtools single-threaded: the threads are Bowtie's (thread counts do not change the alignments)
        'set -o pipefail; bowtie --threads {threads} --sam -x "{params.index}" -q {params.args} "{input.fastq}" 2> "{log}" '
        '| samtools view -bS -o "{output}" -'


rule smrnaseq_mirtop:
    input:
        bam = rules.smrnaseq_align.output,
        hairpin = rules.smrnaseq_hairpin.output,
        gff = rules.mirbase_gff.output
    output:
        gff = join(SMRNASEQ_OUT, 'mirtop', '{sample}', 'mirtop.gff'),
        counts = join(SMRNASEQ_OUT, 'mirtop', '{sample}', 'counts', 'mirtop.tsv')
    params:
        out = join(SMRNASEQ_OUT, 'mirtop', '{sample}'),
        species = SMRNASEQ_SPS
    log:
        join(SMRNASEQ_OUT, 'mirtop', '{sample}', 'mirtop.log')
    container:
        'docker://' + config['docker']['smrnaseq_mirtop']
    shell:
        'mirtop gff --sps {params.species} --hairpin {input.hairpin} --gtf {input.gff} -o {params.out} {input.bam} > {log} 2>&1 && '
        'mirtop counts --add-extra --hairpin {input.hairpin} --gtf {input.gff} --sps {params.species} --gff {output.gff} '
        '-o {params.out}/counts >> {log} 2>&1'


rule smrnaseq_counts:
    input:
        expand(rules.smrnaseq_mirtop.output.counts, sample=SAMPLES)
    output:
        mirna = join(SMRNASEQ_OUT, 'mirna.tsv'),
        joined = join(SMRNASEQ_OUT, 'joined_samples_mirtop.tsv')
    params:
        script = src_gcf('scripts/smrnaseq_counts.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} --tables {input} --joined {output.joined} --mirna {output.mirna}'


rule smrnaseq_all:
    input:
        rules.smrnaseq_counts.output
