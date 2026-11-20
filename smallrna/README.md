# smallrna

**Status: active.** Maintained and covered by a CI dry run and an AnnData regression test. Maintenance status does not mean every kit, route or output listed here has been scientifically validated.

Small RNA sequencing: adapter clipping and trimming, spike-in (calibrator) handling, genome alignment, miRNA, isomiR and tRNA-fragment quantification, and AnnData assembly.
Entry point: [`smallrna.smk`](smallrna.smk). For how BFQ runs a workflow and what it delivers, see the [README](../README.md#how-a-workflow-runs-under-bfq).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: smallrna`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts |
| --- | --- |
| Bioo Scientific NEXTflex Small RNA-Seq Kit v3 | SE |
| Bioo Scientific NEXTflex Small RNA-Seq Kit v4 | SE |
| TruSeq Small RNA Library Prep Kit | SE |

The kit entries set the filter route (`filter.trim.quantifier: cutadapt` with `filter.trim.cutadapt`) and `db.reference_db: ensembl`.
Read 1 reads the RNA in sense orientation (`read_orientation: forward`); no quantifier reads this setting.
NEXTflex v4 offers blockers that deplete abundant RNAs before ligation: a project lists the ones used in `blockers` (`tRNA/YRNA`, `blood miRNA`, `custom`; default none). They are shown in the MultiQC report header and change no processing.

### Filtering with cutadapt (kit default)

Three cutadapt 5.2 passes per sample. [`cutadapt_config.py`](rules/filter/scripts/cutadapt_config.py) gives each rule and the evidence for its lengths and error rates.

1. **3′ adapter** removed.
2. **Filter:** the kit's random ends cut from both ends; reads containing fragments of the 5′ adapter discarded (they pass other filters and are counted as miRNAs); the kit's minimum length, N limit and expected-error limit (`max_ee`).
3. **Artefacts:** reads still containing 3′ adapter or Illumina primer sequence discarded, and reads of the kit's reagent oligos (the TruSeq Stop Oligo).

| Kit | `adapter5` | `five_prime_piece` | `random_ends` | `minimum_length` | `max_n` | `max_ee` | `reagent_oligos` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| NEXTflex v3 | `GTTCAGAGTTCTACAGTCCGACGATC` | 14 | 4 | 15 | 5 | 1 | none |
| NEXTflex v4 | `TCTTTCCCTACACGACGCTCTTCCGATCT` | 15 | 0 | 15 | 5 | 1 | none |
| TruSeq Small RNA | `GTTCAGAGTTCTACAGTCCGACGATC` | 14 | 0 | 15 | 5 | 1 | Stop Oligo (STP) |

All three kits share the 3′ adapter (`adapter`) `TGGAATTCTCGGGTGCCAAGG`. miRge3 and mirflowz take their minimum read length from `minimum_length`.
The passes' reports are merged into `logs/<sample>/<sample>.cutadapt.json` ([`merge_cutadapt.py`](rules/filter/scripts/merge_cutadapt.py)), shown in MultiQC by gcf-multiqc's `gcf_smallrna_cutadapt` module; the per-pass reports are delivered under `provenance/cutadapt/`.

### fastp route

`filter.trim.quantifier: fastp` runs fastp 0.20.1 in two passes, adapter removal and then the project's `filter.trim.fastp.args` ([`cliptrim.smk`](rules/filter/cliptrim.smk) explains why two). No kit selects it; projects initialized before #350 keep it through their own `config.yaml`, which sets both keys. `filter.trim.fastp.args` with another route stops the workflow at load.

The kits' former options:

| Kit | Final options |
| --- | --- |
| NEXTflex v3 | `-f 4 -t 4 -l 16 -q 15 -u 40 -n 5 -e 0 --low_complexity_filter --complexity_threshold 30` |
| NEXTflex v4, TruSeq Small RNA | `-l 15 -q 15 -u 40 -n 5 -e 0 --low_complexity_filter --complexity_threshold 30` |

**Migrating older projects.** Projects with the old sibling key `trim_args` (even empty) or a nonempty `params` stop at load. Remove `trim_args` without appending its `-Q -n 0 -l 0 -G` to `args`: those were first-pass controls, which the workflow now generates. `args` must set a final length policy (`-l N` or `-L`, not both); note that `-Q` also disables N-base filtering. Resumed workdirs keep their captured configuration.

The fastp report `logs/<sample>/<sample>.fastp.json` combines both passes; each pass's native report is under `provenance/fastp/`.

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| gcfdb | `rules/gcfdb.smk` | Reference selection |
| filter | `rules/quant/spikein.smk`, `rules/filter/cutadapt.smk`, `rules/filter/cliptrim.smk`, `rules/filter/spikein_bowtie.smk`, shared `common/` rules | Spike-in sets; FASTQ merging; adapter removal and filtering; spike-in counts and removal |
| align | `rules/align/star.smk` | STAR, the only supported aligner |
| quant | `rules/quant/unitas.smk`, `mirge3.smk`, `mirflowz.smk`, `tesmall.smk`, `tinyrna.smk`, `smrnaseq.smk`, `htseq.smk` | UNITAS (delivered); miRge3, mirflowz, TEsmall, tiny-count, nf-core/smrnaseq, htseq |
| qc | `rules/qc/`, shared `common/` rules | FastQC; miRTrace contamination check and miRNA complexity on the filtered reads |
| bfq | `rules/bfq/` (levels 1–3) | Collect delivery products into `bfq/` |

**`quant.method` and delivery.** `quant.method` selects what the `quant_all` target builds: `unitas`, `mirge3`, `mirflowz`, `tesmall`, `tinyrna`, `smrnaseq`, `htseq`, or several, comma-separated. It does not change delivery: `bfq_all` and `multiqc_report` always deliver the UNITAS products. miRge3, mirflowz, TEsmall, tiny-count, nf-core/smrnaseq and htseq are routes to evaluate (see the [benchmark](benchmark/README.md)).

## Quantification routes

Every route reads the same filtered reads and removes no adapters itself (mirflowz still trims N bases at read ends). There is no read deduplication: `filter.dedup: true` stops the workflow at load. The routes' references are described in [gcfdb](../gcfdb/README.md#small-rna-quantifier-references).

- **UNITAS** 1.8.0 (`gcfntnu/unitas:1.8.0-1`, patched to run offline and fail on errors), the delivered route. Reference dumps exist for human, mouse and rat; another organism stops the workflow at load. With `db.unitas.mirna_db: mirgenedb` miRNAs are annotated against MirGeneDB 3.0 and named by gene and arm (`Let-7-P1b-5p`), so tables are not comparable by name with miRBase runs. Each sample runs in its own work directory ([`run_unitas.sh`](rules/quant/scripts/run_unitas.sh)). The `misc_RNA` class is also split by gene family into `Y_RNA`, `7SK`, `7SL` and `other` ([`unitas_misc_rna.py`](rules/quant/scripts/unitas_misc_rna.py)).
- **miRge3** 0.1.4 (`quant/mirge3/`), with a versioned library (`db.mirge3.release`) and `quant.mirge3.db` (`mirbase` or `mirgenedb`). `miR.Counts.csv` is built by the workflow from miRge3's alignments: exact and isomiR reads together, for miRNAs with at least 2 exact reads in the project ([`mirge3_counts.py`](rules/quant/scripts/mirge3_counts.py)); miRge3's own table is kept as `miR.Counts.mirge3.csv`.
- **mirflowz** 0.10.0 (`quant/mirflowz/`): genome and transcriptome alignment with segemehl, counting miRNAs, isomiRs and pri-miRs from miRBase coordinates. It needs a reference on the miRBase GFF's assembly, so a rat project's mirflowz alignments are on Rnor_6.0. A read aligning equally well to several loci is split between them. Runs nested in `gcfntnu/mirflowz:0.10.0` ([`run_mirflowz.py`](rules/quant/scripts/run_mirflowz.py)), within `quant.mirflowz.mem_gb`.
- **TEsmall** 2.1.0 (`quant/tesmall/`): exact Bowtie 1 alignment to the genome; each read goes to the first class it overlaps by 90% in `quant.tesmall.order` (miRNA first by default; TEsmall's own default puts structural RNA first and takes the reads of miRNAs inside snoRNA genes). `gcfntnu/tesmall:2.1.0-1` with `--no-trim`, so the filtered reads are not trimmed again ([`tesmall.smk`](rules/quant/tesmall.smk)). Its temporary files go to the sample's work directory, and since TEsmall reports success even when a step failed, the job stops on a failed step or an empty annotation (#419). Sample names with dots stop the workflow at load.
- **tiny-count** 1.6.0 (`quant/tinyrna/`), the counter of [tinyRNA](https://github.com/MontgomeryLab/tinyRNA) without the rest of its pipeline: reads collapsed to unique sequences, aligned end to end to the genome with Bowtie 1 (`quant.tinyrna.mismatch`, all alignments of the best stratum) and counted against the miRBase GFF with a Features Sheet written from `quant.tinyrna`: a read counts for a mature miRNA when it starts at one of the 5′ offsets of the annotated 5′ end (default −1, 0, +1) and ends at most `three_prime_extension` beyond its 3′ end (default 4). A sequence's reads are divided between its loci and between the features at a locus (`normalize_by_genomic_hits`, `normalize_by_feature_hits`). Feature IDs are miRBase's, one row per locus (`MIMAT0000416`, `MIMAT0000416_1`); `tinyrna_counts.tsv` joins the samples. `gcfntnu/tiny-count:1.6.0-1`, with a wrapper that checks tiny-count's log and tables, as it exits 0 after an error ([`tinyrna.py`](rules/quant/scripts/tinyrna.py)).
- **nf-core/smrnaseq** 2.4.1's mirtop counting (`quant/smrnaseq/`): reads collapsed to unique sequences, aligned to the organism's miRBase hairpins with Bowtie 1 (up to 2 mismatches in the first 28 bases, none limited beyond; the best alignments, at most 50) and annotated by mirtop; `mirna.tsv` sums each miRNA's isomiRs. The pipeline's steps, images and arguments are reproduced as rules, without Nextflow ([`smrnaseq.smk`](rules/quant/smrnaseq.smk)); on the same reads its tables are identical to the pipeline's. It needs no genome reference.
- **htseq** (`quant/htseq/`), the genome counter of an earlier in-house small-RNA pipeline (#425): reads collapsed to unique sequences with fastx_collapser, those of at least 15 nt aligned end to end to the genome with Bowtie 2 (`-k 10`: up to 10 alignments, mismatches within Bowtie 2's score threshold), and every alignment annotated by htseq-count against the miRBase GFF (union overlap mode, stranded). Only alignments without mismatches count, and a read counts in full for every miRNA name it overlaps at any of its alignments, never divided ([`htseq_counts.py`](rules/quant/scripts/htseq_counts.py)); `htseq_counts.tsv` joins the samples, and each sample also has an isomiR table and statistics. It uses TEsmall's genome and tiny-count's renamed miRBase GFF, without the GFF's rows that have no strand (three human hairpins in miRBase 23.1), which htseq-count's stranded mode cannot place, and its own Bowtie 2 index.

### Spike-in counts

Spike-ins are a project setting: listing sets in `filter.spikein.ref` (gcfdb sets such as `ERCC92` and `smallrna_calibrators`, or a project FASTA) turns them on. By default (`filter.spikein.quantifier: bowtie`) bowtie counts each spike-in from exact, forward-strand matches over the whole read and removes those reads before every other step ([`spikein_bowtie.smk`](rules/filter/spikein_bowtie.smk)); the counts are delivered as `exprs/spikein_counts.tsv` and shown in MultiQC's General Statistics. With `native`, UNITAS and miRge3 estimate the spike-ins themselves from the reads that include them. Switching spike-ins off stops the workflow at load while an earlier run's spike-in outputs remain in `bfq/`.

## Outputs in `bfq/`

`data/tmp/smallrna/bfq/`, as planned by the fixture's dry run (NEXTflex v3). The set is the same for every `quant.method`, because delivery always uses UNITAS.

```text
multiqc_<project>.html, .multiqc_config.yaml
align/<sample>.sorted.bam
exprs/mir_counts.tsv, isomir_counts.tsv, isomir_counts_anno.tsv, trf_counts.tsv, isotrf_counts.tsv, annotations.tsv
exprs/misc_rna_families.tsv  (UNITAS misc_RNA reads by gene family)
exprs/spikein_counts.tsv  (with spike-in sets in filter.spikein.ref)
figs/gene_high_mqc.yaml, pca_mqc.yaml, unitas_misc_rna_mqc.yaml, spikein_mqc.yaml (with spike-ins)
logs/unitas/<sample>  (relative links to quant/unitas/unitas/<sample>)
logs/mirtrace/mirtrace-results.json, mirtrace-report.html, mirtrace-stats-{length,contamination_basic,mirna-complexity}.tsv
logs/<sample>/<sample>.cutadapt.json, _fastqc.html, _fastqc.zip, .Log.final.out  (fastp route: .fastp.json, .fastp.html)
provenance/cutadapt/<sample>.{adapter,filter,artefacts}.{json,log}  (fastp route: provenance/fastp/)
```

- `mir_counts.tsv`, `trf_counts.tsv` and `isomir_counts.tsv` have features as rows and samples as columns (`isomir_counts.tsv` rows are `sequence|miR-name`, annotated in `isomir_counts_anno.tsv`); `isotrf_counts.tsv` and `annotations.tsv` are long tables with a `Sample_ID` column.
- The report's PCA uses log2(CPM + 1) of the miRNAs detected in at least min(3, samples) samples, with a note instead when two components are impossible; `gene_high_mqc.yaml` shows the ten most abundant miRNAs.
- `quant/unitas/adata.h5ad` holds the raw UNITAS miRNA counts with a DESeq2-transformed layer `vsd` (rlog up to 50 samples, otherwise the variance-stabilizing transformation; [`create_anndata.py`](rules/quant/scripts/create_anndata.py)).

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | `unitas`; also `mirge3`, `mirflowz`, `tesmall`, `tinyrna`, `smrnaseq`, `htseq`, comma-separated. Selects the `quant_all` targets, not delivery | [smallrna.config](smallrna.config) |
| `filter.trim.quantifier` | `cutadapt`; also `fastp` or `skip` (reads already trimmed) | Kit entry |
| `filter.trim.cutadapt` | Per kit (table above) | Kit entry |
| `filter.trim.fastp.args` | None; required with the fastp route only | Project |
| `filter.spikein.ref` | `null` (no spike-ins); spike-in sets or project FASTA paths, comma-separated | `smallrna.config` |
| `filter.spikein.quantifier` | `bowtie`; or `native` | `smallrna.config` |
| `filter.contaminants` | Off; `quantifier: unitas` adds the sequences to the UNITAS reference | `smallrna.config` |
| `blockers` | `[]`; NEXTflex v4 only | Project |
| `db.unitas.release`, `db.unitas.mirna_db`, `db.unitas.mirgenedb_release` | `05-11-2026`, `mirbase`, `06-11-2026-mirgenedb` | [main.config](../main.config) |
| `db.mirbase.release` | `23.1`; must equal the UNITAS dump's miRBase release | `main.config` |
| `db.mirge3.release`, `quant.mirge3.db` | `19-11-2026`, `mirbase` | `main.config`, `smallrna.config` |
| `db.mirflowz.ensembl_release`, `db.mirflowz.custom`, `quant.mirflowz.mem_gb` | 116 for human and mouse, 104 for rat; none; 60 | `main.config`, `smallrna.config` |
| `quant.tesmall.order`, `max_length`, `mismatch`, `maxaln` | `miRNA,hairpin,structural_RNA,exon,TE,intron,piRNA_cluster`, 36, 0, 100 | `smallrna.config` |
| `quant.tinyrna.mismatch`, `five_prime_offsets`, `three_prime_extension`, `strand`, `length`, `normalize_by_genomic_hits`, `normalize_by_feature_hits` | 0, `-1,0,1`, 4, `sense`, `16-28`, true, true | `smallrna.config` |
| `align.genome.aligner`, `align.star.args` | `star` (the only accepted value); STAR options for short reads | `smallrna.config` |

Unsupported values of these keys stop the workflow at load, except the tiny-count counting rules (`five_prime_offsets`, `three_prime_extension`, `strand`, `length`), which the `tinyrna_features` job checks when it writes the Features Sheet. Precedence between these sources is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `smallrna/.tests`: synthetic dry-run fixture (NEXTflex v3, `mus_musculus`, four samples). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-small-rna.yml`).
- `smallrna/.tests/test_*.py` and the smallrna references' tests in `gcfdb/.tests/`: each route, the reference checks, the delivery scripts and the load errors, on synthetic data. Tests needing a tool (fastp 0.20.1, cutadapt 5.2, bowtie, samtools) skip without it. CI runs the AnnData regression and the report-figure tests.
- Commands and environments: [development guide](../docs/development.md#run-the-checks).

A dry run shows only that the job graph can be built. See [what the checks prove](../docs/development.md#what-the-checks-prove).

## More documentation

- [`methods.txt`](methods.txt): methods text for reports. It is not read by the workflow.
- [Quantification benchmark](benchmark/README.md)
- [Architecture](../docs/architecture.md), [configuration](../docs/configuration.md), [reference data](../gcfdb/README.md)

## Credits

Written by Arnar Flatberg ([flatberg](https://github.com/flatberg)) at the Genomics Core Facility, [Norwegian University of Science and Technology, NTNU](http://www.ntnu.no). Licence: [LICENSE](LICENSE).
