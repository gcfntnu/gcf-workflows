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

The kit entries set `filter.trim.quantifier: cutadapt` with `filter.trim.cutadapt` (since [#350](https://github.com/gcfntnu/gcf-workflows/issues/350); before, fastp with `filter.trim.fastp.args`) and `db.reference_db: ensembl`.

**Kit metadata** ([#361](https://github.com/gcfntnu/gcf-workflows/issues/361)). Small-RNA libraries are single-end, and read 1 reads the RNA in sense orientation (the 5′ adapter carries the read 1 primer site), so every kit has `read_orientation: forward` and `library_strand: 'second strand'`; no small-RNA quantifier reads them. NEXTflex v4 offers blockers that deplete abundant RNAs before ligation (kit manual: tRNA/YRNA Blockers, in the kit; Blood miRNA Blockers for hsa-miR-451a, hsa-miR-486-5p and hsa-miR-92a-3p; Custom Small RNA Blockers), so those RNAs are deliberately low in the data. Its kit entry has `blockers: []` (none, the default); a project lists the ones used, any of `tRNA/YRNA`, `blood miRNA` and `custom`. The value is shown as *Small RNA blockers* in the MultiQC report header and changes no processing. Unknown values, or blockers with a kit whose entry has no `blockers` field, stop the workflow at load.

### Filtering with cutadapt (kit default)

The kit default (`filter.trim.quantifier: cutadapt`) since #350, using cutadapt 5.2 only ([#350](https://github.com/gcfntnu/gcf-workflows/issues/350); [`cutadapt_config.py`](rules/filter/scripts/cutadapt_config.py)). It removes reads derived from the 5′ adapter, which the fastp stages pass and UNITAS and miRge3 count as miRNAs: in public plasma data (TruSeq, GEO GSE94582) up to 13% of UNITAS's and 48% of miRge3's miRNA reads in a low-yield library, and miR-3168-5p (plus, for miRge3, miR-138-5p) in every sample. Against the fastp route, delivered counts lose these false miRNAs and gain about 1% of reads (cutadapt's 3′ trimming tolerates indels and avoids fastp's occasional one-base over-trim of miRNAs ending in `G`); per-miRNA counts otherwise change by about 1% (#350).

1. **3′ adapter:** `-a <adapter>;min_overlap=5` (cutadapt aligns with indels, so adapters with a deleted base are trimmed).
2. **Filter:** for kits with random ends (`random_ends`, 4 for NEXTflex v3) that many bases are cut from both ends, here rather than in the first pass, because cutadapt cuts before it searches adapters and would otherwise truncate adapter dimers whose adapter starts within the first bases. Reads containing the 5′ adapter are discarded (`--action=none --discard-trimmed`). The fragments start anywhere in the adapter, usually end before its end and sometimes follow a few other bases, while cutadapt accepts a partial 5′ adapter only when it runs to the adapter's end and an internal one only when complete. So `five_prime_adapters.fa` holds the adapter truncated at every end position (non-internal `X` records, matched only at the read start, minimum overlap 12 nt) and every piece of it of at least `five_prime_piece` nt (matched anywhere; 14 for the NEXTflex v3/TruSeq adapter, 15 for NEXTflex v4), each with cutadapt's 10% error rate set explicitly: one error (mismatch or indel) in matches of 12–19 nt, two in matches of 20–26 nt. One error at any length cannot be set for partial matches, and allowing one error, or no indels, kept real artefacts such as adapter fragments with a two-base deletion. The lengths matter: a 10-nt overlap also removed miR-22-5p and miR-139-5p, whose first 11 bases match the adapter with one mismatch, and 12–13-nt pieces anywhere match miR-625-5p; as set, the rule matches only miR-3168-5p among the 2,789 human mature miRNAs of miRBase 23.1. The NEXTflex v4 adapter needs 15-nt pieces: its 14-nt pieces match hsa-miR-6865-3p with one error, while from 15 nt nothing in human mature miRNAs or the small-RNA references matches. 15 nt for the v3/TruSeq adapter would keep artefacts (0.3% of kept reads in a low-yield public plasma library, mostly the adapter's last 11–13 nt after a few other bases, such as `AAAAGTCCGACGATC`). On public NEXTflex v4 plasma (GEO GSE317995, four runs) the v4 settings found the 3′ adapter in 99% of reads and removed 0.02–0.08% of reads as 5′ adapter artefacts; 14-nt pieces would have removed only 3–7 more reads, none a miRNA. Then the kit's minimum length (`-m`), N limit (`--max-n`) and expected-error limit (`--max-ee`, the sum of error probabilities from the base qualities).
3. **3′ artefacts** ([#353](https://github.com/gcfntnu/gcf-workflows/issues/353)): reads still containing 3′ adapter or Illumina primer sequence are discarded in a third pass (`three_prime_artefacts.fa`, `--action=none --discard-trimmed`): mostly inserts whose 3′ adapter has too many sequencing errors for pass 1 (in plasma often the RNY4 Y RNA fragment), adapter dimers with an indel in the adapter, and primer-only reads. With that much adapter left no quantifier annotates them, so no count changes; the report names them. The records are every piece of at least 16 nt of the kit's 3′ adapter followed by the start of Illumina's index-read primer (`…AACTCCAGTCAC`), and of P7 (`ATCTCGTATGCCGTCTTCTGCTTG`), with the 10% error rate; the Illumina sequences are workflow constants, the same for every kit. With 14-nt pieces, pieces at the adapter's end matched human 18S rRNA, hsa-miR-2116-3p and a misc_RNA; from 16 nt nothing in human mature miRNAs or Ensembl/GtRNAdb/SILVA small-RNA references matches. In public plasma data the pass removed 0.07–0.92% of reads (TruSeq) and 0.3–0.65% (NEXTflex v4, GSE317995); on a NEXTflex v3 pool run 0.18%, with no reads exactly matching pool sequences lost. A few hundred reads per library keep a short, heavily damaged adapter end (e.g. a 4-base deletion in 19 nt) that neither pass 1 nor these records match. It is a separate pass so that its count is exact: cutadapt matches adapters before it applies the length and N filters.

Settings per kit in `libprep.config` under `filter.trim.cutadapt`: `adapter5` (the 5′ adapter), `five_prime_piece` (at least 14), `random_ends`, `minimum_length`, `max_n`, `max_ee` (`null` for none); the 3′ adapter is the kit's `adapter`. Every small-RNA kit has them:

| Kit | `adapter5` | `five_prime_piece` | `random_ends` | `minimum_length` | `max_n` | `max_ee` |
| --- | --- | --- | --- | --- | --- | --- |
| NEXTflex v3 SE | `GTTCAGAGTTCTACAGTCCGACGATC` | 14 | 4 | 15 | 5 | 1 |
| NEXTflex v4 SE | `TCTTTCCCTACACGACGCTCTTCCGATCT` (NEXTflex 5′ Adapter v4, the Illumina Read 1 sequence; kit manual, Appendix A) | 15 | 0 (v4 has no random bases) | 15 | 5 | 1 |
| TruSeq Small RNA SE | `GTTCAGAGTTCTACAGTCCGACGATC` | 14 | 0 | 15 | 5 | 1 |

All three kits share the 3′ adapter `TGGAATTCTCGGGTGCCAAGG`. The minimum length is 15 nt for every kit (NEXTflex v3 had 16 nt with fastp, #215). **Quality:** `max_ee: 1` discards reads expected to hold one or more sequencing errors, a read-level filter like the fastp route's `-q 15 -u 40`. On 11 public runs (TruSeq and NEXTflex v3 pools, TruSeq plasma; HiSeq and NextSeq) it removed 0–1.3% of reads (mostly one-off low-quality reads, some with N bases or polyG), lost no read exactly matching a human mature miRNA in plasma and at most 2 exact pool reads per run, and changed no miRNA count; on four public NEXTflex v4 plasma runs it removed 0.03–0.08% of reads and no exact miRNA read. The manufacturer's v4 trimming instructions use `--quality-cutoff 20` (3′ end trimming) and a 16-nt minimum instead; on the same runs `-q 20` lost 5–40 exact pool reads per run and gained nothing (#350). A settings table without `filter.trim.cutadapt` stops the workflow at load. miRge3 and mirflowz take their minimum length from `minimum_length`.

**fastp settings.** The kits no longer set `filter.trim.fastp.args`. Projects initialized before #350 carry the kit's previous route and fastp options in their `config.yaml`, which wins over the kit entry, so they keep running fastp. `filter.trim.fastp.args` with any other route (`cutadapt` or `skip`) stops the workflow at load instead of being ignored: remove it, or set `filter.trim.quantifier: fastp`.

**Reports.** Each pass writes a cutadapt JSON report (`filter/cutadapt_stages/<sample>.{adapter,filter,artefacts}.json`, with the console output in `.log`). [`merge_cutadapt.py`](rules/filter/scripts/merge_cutadapt.py) merges them into one report in cutadapt's schema, delivered as `logs/<sample>/<sample>.cutadapt.json` instead of the fastp report: input reads from the first pass, filtered categories of the first two passes summed (5′ adapter artefacts as `discard_trimmed`), the third pass's discarded reads as `gcf_three_prime_artefacts`, output reads from the third pass, trimmed-length statistics of the 3′ adapter, and a record of the three passes (`gcf_passes`). The pass reports are delivered as `provenance/cutadapt/`. The MultiQC report reads the merged reports with gcf-multiqc's `gcf_smallrna_cutadapt` module (gcfntnu/gcf-multiqc#9; selected in the smallrna module list, #358), in fastp's place: every filtering category named, `discard_trimmed` as *5′ adapter artefacts*, and General Statistics with the reads before filtering, % passing filter and % 3′ adapter, placed as fastp's. gcf-bfq's analysis QC reads the merged report (gcfntnu/gcf-bfq#187).

### fastp route (`filter.trim.quantifier: fastp`)

Optional since #350: no kit sets it, so a project that selects fastp also sets `filter.trim.fastp.args` (projects initialized before #350 have both). The route runs **adapter removal, then final trimming/filtering**, with fastp 0.20.1.
The workflow generates the adapter-only options `-Q -L -G`: no quality/N-base filtering, optional length filtering or polyG trimming in that pass.
Empty/null reads are still rejected by fastp and counted in that stage's `too_short_reads`; `-L` does not preserve empty reads.
Only the second invocation receives `filter.trim.fastp.args`; the workflow supplies `-A` there to prevent repeated adapter removal.
There is no configurable first-stage option string, and no small-RNA use of the shared rule's `params` key.

The options the kit entries set before #350:

| Kit | Fixed clipping after adapter removal | Minimum final insert |
| --- | --- | --- |
| NEXTflex v3 SE | Four bases from each exposed end (`-f 4 -t 4`) | 16 nt |
| NEXTflex v4 SE | None | 15 nt (existing policy preserved) |
| TruSeq Small RNA SE | None | 15 nt (existing policy preserved) |

Every former preset sets the same final quality/N-base and low-complexity thresholds:
`-q 15 -u 40 -n 5 -e 0 --low_complexity_filter --complexity_threshold 30`.
This rejects reads with more than 40% of bases below Q15 or more than five N bases, with no average-quality requirement;
low complexity is filtered at fastp's 30% adjacent-base-change threshold.
Final polyG handling retains fastp 0.20.1's instrument-name detection (enabled automatically for NextSeq/NovaSeq), after adapter removal.
`-G` in a deliberate final override disables that behavior. Presets do not enable quality end cutting or polyX trimming.

**Why retain two passes for v4 and TruSeq?** A single invocation does not preserve the existing polyG operation order.
fastp trims polyG before adapters within an invocation; a polyG tract exposed only after adapter removal therefore behaves differently.
The synthetic regression uses a NextSeq-style header and demonstrates a retained-read difference between the proposed one-pass route and the existing route.
For v3, fixed end trimming also precedes adapter removal within an invocation, so concatenating its options would clip the wrong 3′ bases.
See [fastp 0.20.1 processing order](https://github.com/OpenGene/fastp/blob/v0.20.1/src/seprocessor.cpp).
This bounded correction keeps both stages for all kits; a one-pass change needs a separately validated scientific policy.

The v3 16 nt minimum was the workflow requirement in [#215](https://github.com/gcfntnu/gcf-workflows/issues/215); the cutadapt route uses 15 nt for every kit.

### Existing project overrides

Projects initialized from #215 until #350 have the kit's final options under the string key `filter.trim.fastp.args`.
Projects retaining the old sibling `trim_args` (even empty), or a nonempty `params`, fail early with migration guidance instead of silently changing the stages.
For an intentional migration, inspect the old effective configuration and remove `trim_args`; do **not** append its `-Q -n 0 -l 0 -G` to the final options.
Those were first-pass controls. The workflow now supplies equivalent adapter-only controls independently of final filtering.
Set `args` to the applicable preset above, including v3's `-f 4 -t 4 -l 16`, or preserve explicitly reviewed project-specific final overrides.
The argument string must specify a final length policy (`-l N` or `-L`, not both); use `-l 0` or `-L` only for deliberately disabled minimum-length filtering.
In fastp, `-Q` also disables N-base filtering, regardless of `-n`; document that choice when overriding the final policy.
Resuming an older workdir uses its captured workflow/configuration; deploying this revision does not migrate it automatically.

### fastp reports

With the fastp route, the `<sample>.fastp.json` delivery path now contains reconciled original-input/final-output metrics, both commands in execution order,
and a `processing_stages` array with each stage's summary and rejection counts.
The merge rejects inconsistent stage accounting or a broken read/base-count chain.
Raw-input duplication and adapter statistics come from the first pass; final retained-read statistics come from the second.
The existing `<sample>.fastp.html` is a combined stage/count summary rather than the second pass's native HTML.
Native reports remain available as `provenance/fastp/<sample>.adapter.{json,html}` and `<sample>.filter.{json,html}`.
They are outside MultiQC's `logs/` scan so the report sees one combined sample, without duplicate stage samples.

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| gcfdb | `rules/gcfdb.smk` | Reference selection |
| filter | `rules/quant/spikein.smk`, `rules/filter/cutadapt.smk`, `rules/filter/cliptrim.smk`, `rules/filter/spikein_bowtie.smk`, shared `common/` rules | Spike-in sets; FASTQ merging; adapter removal and filtering (cutadapt by default, fastp optional); spike-in counts per spike-in and optional removal (bowtie) |
| align | `rules/align/star.smk` | STAR (smallrna options); the only supported genome aligner |
| quant | `rules/quant/unitas.smk`, `mirge3.smk`, `mirflowz.smk`, `tesmall.smk` | UNITAS (including its variance-transformation rules); miRge3; mirflowz; TEsmall |
| qc | `rules/qc/`, shared `common/` rules | FastQC (shared); miRTrace on the filtered reads the quantifiers read, without adapter trimming (no adapter, `-p illumina`): its contamination check and miRNA complexity are shown in the MultiQC report, its QC plot, RNA categories and read lengths are not (they repeat the filter, UNITAS and FastQC sections). miRTrace's fixed filters (under 18 nt, low complexity, low PHRED) only leave reads out of its own statistics. Delivered to `bfq/logs/mirtrace/` ([#352](https://github.com/gcfntnu/gcf-workflows/issues/352)) |
| bfq | `rules/bfq/` (levels 1–3) | Collect delivery products into `bfq/` |

**`quant.method` and delivery.** `quant.method` selects what the separate `quant_all` target builds: `unitas`, `mirge3`, `mirflowz`, `tesmall`, or several of them (comma-separated). Any other value stops `quant_all` with an error. It does **not** change delivery: `bfq_all` and `multiqc_report` always deliver the UNITAS products (`rules/bfq/bfq_level2.smk`), whatever `quant.method` says.

**Quantifier status.** UNITAS is no longer distributed upstream, and its download server is gone. The workflow uses `gcfntnu/unitas:1.8.0-1`, UNITAS 1.8.0 patched to run offline and to exit non-zero on errors, with reference dumps built from Ensembl, GtRNAdb and miRBase (`db.unitas.release`). Dumps exist only for `homo_sapiens`, `mus_musculus` and `rattus_norvegicus`; any other organism stops the workflow at load. UNITAS remains the delivery route until a maintained quantifier replaces it. miRge3 0.1.4 is kept as an alternative `quant_all` route, with versioned libraries (`db.mirge3.release`, see [miRge3 libraries](#mirge3-libraries)), and mirflowz 0.10.0 and TEsmall 2.1.0 are added as ones to evaluate ([mirflowz](#mirflowz), [TEsmall](#tesmall)). miRge (v2), QuickMIRSeq and SeqBuster were removed.

**UNITAS execution.** Each sample runs in its own work directory, `quant/unitas/work/<sample>/`, with the reference dump linked in ([`run_unitas.sh`](rules/quant/scripts/run_unitas.sh)); `seqmap.exe` comes from the image. Exactly one UNITAS result directory must appear there; it is moved to `quant/unitas/unitas/<sample>/` and the work directory is removed. Concurrent or similarly named samples and leftovers from interrupted runs therefore cannot mix. UNITAS console output goes to `quant/unitas/logs/<sample>.log`. The job fails, and prints the end of that log, when UNITAS exits non-zero, reports `ERROR!` or writes no `unitas.annotation_summary.txt`. Spike-in and contaminant sequences are passed to UNITAS with `-refseq` only when `filter.spikein.quantifier` or `filter.contaminants.quantifier` is `unitas`.

## Outputs in `bfq/`

`data/tmp/smallrna/bfq/`, as planned by the fixture's dry run (NEXTflex v3). The set is the same for every `quant.method`, because delivery always uses UNITAS.

```text
multiqc_<project>.html, .multiqc_config.yaml
align/<sample>.sorted.bam
exprs/mir_counts.tsv, isomir_counts.tsv, isomir_counts_anno.tsv, trf_counts.tsv, isotrf_counts.tsv, annotations.tsv
exprs/spikein_counts.tsv  (with spike-in sets in filter.spikein.ref)
figs/gene_high_mqc.yaml, pca_mqc.yaml (a note instead of the PCA when two components are impossible)
figs/spikein_mqc.yaml  (with spike-in sets in filter.spikein.ref: spike-in reads in General Statistics)
logs/unitas/<sample>  (relative links to quant/unitas/unitas/<sample>)
provenance/cutadapt/<sample>.{adapter,filter,artefacts}.{json,log}  (fastp route: provenance/fastp/<sample>.{adapter,filter}.{json,html})
logs/mirtrace/mirtrace-results.json, mirtrace-report.html, mirtrace-stats-{length,contamination_basic,mirna-complexity}.tsv
logs/<sample>/<sample>.cutadapt.json, _fastqc.html, _fastqc.zip, .Log.final.out  (fastp route: .fastp.json, .fastp.html instead of .cutadapt.json)
```

**MirGeneDB dumps.** With `db.unitas.mirna_db: mirgenedb`, UNITAS annotates miRNAs against MirGeneDB 3.0 instead of miRBase: the curated set (human 514 precursors, 559 mature and 468 star sequences, against about 1,900 human miRBase hairpins), with the other 113 MirGeneDB species for cross-species hits. The rest of the dump (Ensembl, GtRNAdb, piRNA clusters, tRF-1, tRNA leaders, SILVA) is the same as in the miRBase dump. miRNAs are named by MirGeneDB gene and arm, without the species code: `Let-7-P1b-5p` and `Let-7-P1b-3p` for MirGeneDB's `Hsa-Let-7-P1b_5p` and `Hsa-Let-7-P1b_3p*`. That is the name UNITAS also gives reads it places on a hairpin arm (isomiRs), so each arm is counted under one name. Tables are not comparable by name with miRBase runs. `mirgenedb_genes` annotates the AnnData from the dump's MirGeneDB GFF: `star` marks the star arm, and `alias` holds the miRBase accession where there is one. UNITAS' `db_versions.info` still calls the database miRBase; the UNITAS log records `UNITAS MirGeneDB,3.0`. miRge3 keeps using `quant.mirge3.db`.

**UNITAS reference checks.** Two rules identify the selected dump and keep its counts joined to a matching annotation; neither downloads or rebuilds a dump already in `EXT_DIR`, and `unitas_db` itself is unchanged.

- `unitas_db_verify` (once per dump, beside it in `EXT_DIR/unitas/<release>/verify/UNITAS_refdump_<organism>.tsv`): when the release publishes `SHA256SUMS` (releases from `04-11-2026` on), the archive is downloaded, its SHA-256 must equal the published one, and the cached dump must hold exactly the archive's files with the same contents. A difference stops the workflow and names the files; the cache is left as it is. Releases without `SHA256SUMS` are recorded as `not verified`, and remain usable.
- `unitas_reference` (per project, `quant/unitas/reference.tsv`): records the dump's miRNA database and release (from `gcf_mirna_db.tsv` where the dump has it, otherwise UNITAS' `db_versions.info`), its Ensembl, GtRNAdb, piRNA cluster, tRF-1, tRNA leader and SILVA versions, whether it has source checksums (`gcf_sources.tsv`), the verification status, and the annotation's release, each with the file it comes from. Records a dump does not have are `NA`; nothing is filled in from current settings. It is an input of `unitas_mirna_anndata` and `unitas_all`.

UNITAS miRNA counts are joined to the feature annotation by name, and names change between releases. On synthetic reads of every human mature miRNA, the `05-11-2026` dump (miRBase 23.1) leaves 6 of 2,730 counted miRNAs unannotated with the miRBase 23.1 GFF (mature sequences without a genome locus there), but 963 (35%) with the miRBase 22 GFF; MirGeneDB counts match none of the miRBase names. `unitas_reference` therefore stops when the dump's miRNA database is not `db.unitas.mirna_db`, when a miRBase dump's release differs from `db.mirbase.release` (the error names both, and suggests `db.mirbase.release` only for releases gcfdb can download), when a MirGeneDB dump disagrees with its own GFF, or when the dump's miRNA release cannot be read. The dump's Ensembl release is recorded, not checked: UNITAS' reference need not match the genome alignment reference.

## Spike-in counts (bowtie)

Spike-in sets in `filter.spikein.ref`, with `filter.spikein.quantifier: bowtie`, the default ([#336](https://github.com/gcfntnu/gcf-workflows/issues/336); [`spikein_bowtie.smk`](rules/filter/spikein_bowtie.smk)). The configured sets become one bowtie reference per project (`filter/spikein/spikein.fa`, DNA, names up to the first space; a name in two sets stops the run). Each sample's filtered reads are matched with bowtie 1.2.3 exactly over their whole length and on the forward strand only (`-v 0 --norc`), so a read counts for a spike-in only when it is an exact piece of it; a read matching two spike-ins counts as `__ambiguous` (bowtie reports every alignment, `-a`: `-k` caps alignment positions, not distinct spike-ins, so a read matching one spike-in twice and another once would be counted for the first). Every sample gets a table with every spike-in, zeros included (`filter/spikein/bowtie/<sample>.tsv`), and the project a matrix of them, delivered as `exprs/spikein_counts.tsv`, so tables of different samples and projects line up. The MultiQC General Statistics table shows each sample's spike-in reads and their percentage of the filtered reads (`figs/spikein_mqc.yaml`, [`spikein_mqc.py`](rules/bfq/scripts/spikein_mqc.py)), placed after the filtering columns; with `quantifier: native` the same columns come from UNITAS's `Spikein` class, which allows mismatches and 3′ tails (on synthetic reads it counted reads with one mismatch that bowtie does not). The route is a prefilter, and runs whenever spike-in sets are listed, with either quantifier: a step sees spike-in reads only when it estimates them itself. Every step without spike-in support (mirflowz, miRTrace, STAR, FastQC, FastQ Screen) reads the reads that matched no spike-in (bowtie's `--un`; `get_filtered_fastq`). With `bowtie`, UNITAS and miRge3 read them too and estimate no spike-ins (UNITAS gets no spike-in reference, miRge3 no spike-in index). With `native`, UNITAS and miRge3 read the reads with the spike-ins and estimate them (`-refseq`, spike-in index), and bowtie's table is a cross-check. The route replaces the former `rules/filter/calibrators.smk`, which the entry point did not load and which used other set names and index paths and seed-only matching.

## miRge3 libraries

`db.mirge3.release` selects the miRge3 library; `quant.mirge3.db` (`mirbase` or `mirgenedb`) selects the miRNA database within it. The default release, `09-11-2026`, holds libraries for human, mouse and rat built by [gcfntnu/gcf-docker#41](https://github.com/gcfntnu/gcf-docker/pull/41) from miRBase 23.1, MirGeneDB 3.0 and Ensembl 116, in the format of the SourceForge libraries. Each library records its sources: `gcf_mirna_db.tsv` (miRNA releases), `gcf_manifest.tsv` (releases, assemblies, tools) and `gcf_sources.tsv` (every download with URL and SHA-256).

- `mirge3_library` (once per release and organism, [`gcfdb/scripts/mirge3_library.py`](../gcfdb/scripts/mirge3_library.py)): downloads `miRge3_lib_<organism>.tar.gz` from `<winecellar.url>/mirge3/<release>/`. The archive's SHA-256 must equal the release's `SHA256SUMS` before it is unpacked, and the library must hold every file miRge3 opens for both databases (and the tRNA fragment annotation for human), with all six files of every Bowtie index, plus the three provenance files. It is then installed as `EXT_DIR/mirge3/<release>/<org>/` and recorded in `EXT_DIR/mirge3/<release>/verify/miRge3_lib_<organism>.tsv`. A missing `SHA256SUMS`, a checksum or file mismatch, or an existing `<org>/` without that record stops the workflow; each release has its own directory, so other caches are not touched.
- `db.mirge3.release: sourceforge` selects the unversioned SourceForge library at its old place, `EXT_DIR/mirge3/<organism>/<org>/`, for retained analyses. Existing copies are used as they are; the library has no checksums or provenance.
- `mirge3_reference` (per project, `quant/mirge3/reference.tsv`): stops when the library lacks the files of `quant.mirge3.db` (including any of the six files of its Bowtie indexes), or when its `gcf_mirna_db.tsv` does not list that database. It records the database and its release, the Ensembl release and assemblies from the manifest, the source checksums and the archive check, each with the file it comes from; the SourceForge library's are `NA`. miRge3 counts are not joined to another annotation, so no other release is checked. It is an input of `mirge3_quant` and `quant_all`.

**Input.** miRge3 reads the shared filtered reads, as UNITAS and mirflowz do: the kit's adapter removal, fixed clipping, length policy and quality and complexity filtering are done before it (with `filter.trim.quantifier: skip`, the merged reads unchanged: then no step removes adapters for any quantifier, so `skip` is only for reads trimmed before they reach the workflow). So miRge3 gets no `--adapter` (cutadapt then trims nothing), no `--nextseq-trim` (a polyG tail follows the 3′ adapter and is removed with it; the fastp route also trims polyG) and no `--uniq-mol-ids`, and its minimum read length is the kit's final minimum (`filter.trim.cutadapt.minimum_length`; with the fastp route `-l` of `filter.trim.fastp.args`, 1 with `-L`, where the filter keeps every non-empty read) instead of its own 16 nt. There is no deduplication: the NEXTflex v3 randomized 4N ends reduce ligation bias but, used as molecule identifiers (8 random nucleotides), saturate for abundant miRNAs; the filter clips them. Until #342, miRge3 read the merged raw reads, trimmed them itself and collapsed NEXTflex v3 reads by their 4N ends, so its counts for existing projects change; delivery (UNITAS) does not.

**Counting.** miRge3 counts a miRNA's exact and isomiR reads together. It has an option hidden from `--help`, `--crThreshold` (default 0.1), that sets a miRNA's count in a sample to 0 when its exact reads are at most that share of its isomiR reads; the workflow passes `--crThreshold 0`, so isomiR-rich miRNAs keep their reads (#347). Its table also has a fixed rule: a miRNA with fewer than 2 exact reads in a sample is 0 there. The workflow therefore builds `miR.Counts.csv` itself from miRge3's `mapped.csv` (`mirge3_counts`, [`mirge3_counts.py`](rules/quant/scripts/mirge3_counts.py), #349): each miRNA's exact and isomiR reads per sample, with names folded by the library's merge file as miRge3 does, for miRNAs with at least 2 exact reads in the project (all samples together); other miRNAs are 0. Requiring exact reads keeps miRge3's protection against isomiR-only assignments, which on the synthetic equimolar pool (Giraldez et al. 2018) were mostly reads of non-human pool sequences: counting every assigned read detected about 220 miRNAs per sample that are not in the pool, against 10–36 with an exact-read rule. Applying it per project instead of per sample recovered 26–28 more pool miRNAs per sample than miRge3's table. miRge3's own table is kept as `miR.Counts.mirge3.csv`. UNITAS has no such rules.

**Spike-ins.** Spike-ins are a project setting, not a library prep kit setting: the lab adds them to a project's libraries independently of the kit, so no kit in `libprep.config` sets them, and a project that uses them lists them in `filter.spikein.ref` in its own config; that one key turns them on (default `null`: no spike-in step runs), and `filter.spikein.quantifier` is needed only for `native`. Switching spike-ins off in a project that had them stops the workflow at load while `bfq/figs/spikein_mqc.yaml` or `bfq/exprs/spikein_counts.tsv` from the earlier run remain, since the report and delivery would still include them; delete them first. `filter.spikein.ref` lists the sets, gcfdb sets by name (sequences in `EXT_DIR/spikein/fasta/<set>.fa`) or a project's own FASTA by its path (`.fa`/`.fasta`, relative to the working directory, checked at load: named records with A/C/G/T/U sequences; used as a copy in DNA letters, `filter/spikein/sets/custom<n>.fa`, so RNA-letter files work for every route), and every route uses the same sequences ([`spikein.smk`](rules/quant/spikein.smk)). miRge3 0.1.4 counts spike-ins (`--spikeIn`) only from the library's `index.Libs/<org>_spike-in` index and takes no sequences itself. The published library's spike-in index is therefore removed when it is installed, and the shared library is not changed after that. Each project gets its own library, `quant/mirge3/lib/<org>/`, made by `mirge3_project_library` ([`mirge3_project_library.py`](rules/quant/scripts/mirge3_project_library.py)): links to every file of the shared library, plus a spike-in index built from the project's sets, or none when no sets are listed; `--spikeIn` is passed only then, and miRge3's `--libraries-path` is `quant/mirge3/lib`. Projects therefore cannot use each other's spike-ins, and a change of sets rebuilds the project library. `quant/mirge3/spikein.tsv` records the shared library, the sets, their sequence count and the SHA-256 of the combined FASTA. miRge3 aligns spike-ins last, exactly, and only for reads nothing else matched, and reports one `Spike-in` total per sample in `annotation.report.csv`; counts per spike-in come from the bowtie route or UNITAS (`filter.spikein.quantifier: bowtie` or `native`). miRge3 gets the spike-in index only with `native`.

**Changes from the SourceForge libraries.** The new libraries are a scientific reference update; counts are not comparable by name with SourceForge runs, and unchanged counts are not expected. From gcf-docker#41's comparison:

- miRBase 22 to 23.1: renames (`hsa-miR-1-3p` to `hsa-miR-1-2-3p`, `hsa-miR-107` to `hsa-miR-107-3p`) and mature boundaries moved by 1-2 nt (`mmu-let-7g-5p` gains a 3' T). Same name with changed sequence / renamed / dropped / added: human 451 / 743 / 120 / 252, mouse 373 / 269 / 92 / 93, rat 220 / 117 / 30 / 132.
- Merge groups are miRNAs with identical mature sequences and a miRNA's SNP variants (for `--AtoI`, from dbSNP). The SourceForge libraries also merged similar but not identical sequences by an undocumented rule; such miRNAs (e.g. `miR-320a/b/c/d/e`) are now reported separately. Human miRBase merge groups: 284 to 226.
- MirGeneDB 2.0 to 3.0: renamed families (`Hsa-Let-7-P2_5p` to `Hsa-Let-7-P1b_5p`), and the `-v1`/`-v2` arm variants are gone (human mature 1171 to 1027).
- Other RNA from Ensembl 116 (human mRNA 170,448 to 465,253 sequences; `rRNA_pseudogene` counted as rRNA, as in the old library). Mouse and rat tRNAs from GtRNAdb; the human tRNA set and its tRNA fragment files are unchanged. The old rat library was far older (Rnor_5.0 transcripts).
- No genome pickles, so `--novel-miRNA` is not available (not used by the workflow).
- miRge3 0.1.4 writes `miRBase22` or `MirGeneDB2.0` into its GFF output whatever the library; `reference.tsv` records the actual release.

**Delivered UNITAS tables.** `mir_counts.tsv`, `trf_counts.tsv` and `isomir_counts.tsv` have features as rows and samples as columns; `isomir_counts.tsv` rows are `sequence|miR-name`, annotated in `isomir_counts_anno.tsv`. `isotrf_counts.tsv` (tab-separated) and `annotations.tsv` are long tables with a `Sample_ID` column. Every `annotations.tsv` row names its parent levels (`Level1`, `Level2`, `Level3`). All tables are sorted, so repeated runs give identical files.

## mirflowz

`quant.method: mirflowz` (alone or combined, e.g. `unitas,mirflowz`) runs [mirflowz](https://github.com/zavolanlab/mirflowz) 0.10.0 as an alternative `quant_all` route, to evaluate it alongside UNITAS and miRge3; delivery does not use it. mirflowz maps the filtered reads to the genome and transcriptome (segemehl, oligomap) and counts miRNAs, isomiRs and pri-miRs from the alignments' overlap with miRBase coordinates, so unlike UNITAS and miRge3 it depends on the genome assembly. It runs unmodified in `gcfntnu/mirflowz:0.10.0` (gcfntnu/gcf-docker#42, `sha256:5841b820d2ae0f008ca6e0b668069b68eacad32756040c84d9b379c435085e45`), which holds every tool it calls, so the nested run needs no inner containers.

**Reference** ([`gcfdb/mirflowz.db`](../gcfdb/mirflowz.db)), once per reference in `EXT_DIR/mirflowz/mirbase-<release>_ensembl-<release>/<organism>/`:

- Genome, GTF and miRBase GFF must be on one assembly. The Ensembl release is fixed per organism in `db.mirflowz.ensembl_release`, on the assembly of the miRBase 23.1 GFF: human GRCh38 and mouse GRCm39 (Ensembl 116), rat Rnor_6.0 (Ensembl 104, the last release on it). It is independent of `db.ensembl.release`: a rat project keeps its own genome alignment on GRCr8, while mirflowz's BAM files are on Rnor_6.0. A release on another assembly stops the workflow at load; the downloaded GFF's `genome-build-id` is checked again when the reference is built. There is no automatic choice and no liftover.
- The project's Ensembl genome and GTF are reused when `db.reference_db` is `ensembl` and its release and assembly are mirflowz's; otherwise they are downloaded into `sources/`. The GFF comes from gcfdb's miRBase getter.
- The chromosome map from the GFF's names to the genome's (`chr1` to `1`, `chrM` to `MT`, `chrUn_JH584304` to `JH584304.1`; mouse GFFs already use Ensembl names) is built from the files and stops on any GFF chromosome it cannot place ([`mirflowz_chr_map.py`](../gcfdb/scripts/mirflowz_chr_map.py)). The rat GFF's `901866519`, the NCBI GI of a *Plasmodium vivax* contig, not a rat sequence, is removed explicitly with its three features (`rno-mir-16-2` and its arms).
- mirflowz's own preparation rules build the segemehl indexes and processed genome, transcriptome and annotations into `prepared/` ([`mirflowz_prepare.py`](../gcfdb/scripts/mirflowz_prepare.py)), recorded in `mirflowz_reference.tsv` with releases, assembly and input checksums. mirflowz's top-level Snakefile deletes this directory after every run, so the script runs its rule files under its own small Snakefile instead. The GFFs miRBase now serves repeat a mature miRNA's `ID` at each of its loci, which mirflowz's annotation step rejects; the GFF given to mirflowz restores miRBase's older layout, later occurrences suffixed `_1`, `_2`, ... (`Alias` unchanged), and the record counts them.
- `db.mirflowz.custom` (`name`, `genome`, `gtf`, `mirna_gff`, all required together) replaces all three inputs, e.g. with a reduced reference; it is stored as `EXT_DIR/mirflowz/custom-<name>/<organism>/`, and no assembly is checked.

**Run** (`mirflowz_quant`, [`run_mirflowz.py`](rules/quant/scripts/run_mirflowz.py)): mirflowz's mapping and quantification rules run in a work directory that links the prepared reference. mirflowz always quality-filters and adapter-trims its input; on the shared filtered FASTQ this is made a no-op where mirflowz allows it: the reads are passed as FASTA, which skips its quality filter (it filters only FASTQ), `overlap` equal to the adapter length (only a complete adapter is trimmed, so a read whose real 3' end matches the start of the adapter keeps it; mirflowz's default of 3 would cut such ends), `minimum_length` the kit's final minimum (`-l` of `filter.trim.fastp.args`, or `filter.trim.cutadapt.minimum_length`), and no limit on N bases. mirflowz still trims N bases at read ends (`cutadapt --trim-n`, not configurable). mirflowz reads its sample table with pandas type inference, so sample names such as `1` or `01` would turn into numbers and collide; samples are passed to it as `gcf_sample_<n>` and get their own names back in the count tables and alignment directories (`mirflowz_samples.tsv` keeps both).

**Resources.** mirflowz declares threads and memory as custom resources (e.g. `threads=8, mem=50` GB for genome mapping), which Snakemake's core count does not limit. Both nested runs cap every rule's threads at the outer rule's threads and its memory at the outer job's effective `mem_mb` in whole GB (requested as `quant.mirflowz.mem_gb`, default 60, and lowered by Snakemake's `--resources` or `--set-resources` when those apply), and limit the totals to the same budget, so concurrent jobs cannot together exceed what the outer job was given.

mirflowz divides the count of a read that aligns equally well to several loci between them, so a sequence repeated in the genome contributes only part of its reads to an annotated miRNA (UNITAS and miRge3 match by sequence instead).

**Outputs** in `quant/mirflowz/`: `all_mirna_counts.tab` (miRNAs and isomiRs), `all_pri-mir_counts.tab`, `alignments/<sample>/` (mirflowz's per-sample SAM and BAM files) and the generated `mirflowz_config.yaml` and `mirflowz_samples.tsv`.

Full-genome references have not been built yet: segemehl needs about 60 GB of memory for a human or mouse genome index, so building them is a separate task on a large machine. On a reduced reference (mouse chromosome 19, 61 Mb), preparation took 1.5 minutes with 1.2 GB of memory, and quantification of four small samples 20 seconds with 1 GB.

## TEsmall

`quant.method: tesmall` (alone or combined, e.g. `unitas,tesmall`) runs [TEsmall](https://github.com/mhammell-laboratory/TEsmall) 2.1.0 as a `quant_all` route ([#374](https://github.com/gcfntnu/gcf-workflows/issues/374)); delivery does not use it. TEsmall aligns reads to the genome with Bowtie 1 (exact by default, at most 100 alignments, after removing reads that match the rDNA) and assigns each read to the first class it overlaps by at least 90% of the read, in the order structural RNA, miRNA, hairpin, exon, TE, intron, piRNA cluster. Reads hitting several features within a class are split equally between them, and TEsmall rounds the per-feature counts, so shares below 0.5 are lost. Reads ending in CCA first go through its 3′ tRF step against the tRNA genes. A 3′ isomiR that extends past the mature end is counted under the precursor (`hairpin`), not the mature miRNA.

**Database** ([`gcfdb/tesmall.db`](../gcfdb/tesmall.db), [`tesmall_db.py`](../gcfdb/scripts/tesmall_db.py)), once per reference in `EXT_DIR/tesmall/<mirflowz reference name>/<organism>/`, on the genome of the [mirflowz](#mirflowz) reference (the same Ensembl release and assembly, chromosome map and miRBase release, or `db.mirflowz.custom`), instead of TEsmall's prebuilt databases. `tesmall_reference.tsv` records the sources, checksums and counts per class.

- `miRNA`, `hairpin`: miRBase GFF mature and precursor features, by miRBase name.
- `structural_RNA`: Ensembl genes of the biotypes rRNA, Mt_rRNA, Mt_tRNA, snRNA, snoRNA, scaRNA, misc_RNA, vault_RNA, ribozyme, sRNA and scRNA, and GtRNAdb tRNA genes (GRCh38, GRCm39, Rnor_6.0).
- `exon`, `intron`: exons and the gaps between them of all other Ensembl genes except miRNA genes.
- `TE`: UCSC RepeatMasker elements of the classes LINE, SINE, LTR, DNA, RC and Retroposon.
- `piRNA_cluster`: the piRNAclusterDB clusters of the configured UNITAS reference dump (`piRclDB`), placed on the genome by unique exact matches of anchors from both ends of each cluster sequence, so the mouse clusters, published on GRCm38, land on GRCm39 without a liftover. Clusters that cannot be placed uniquely and whole are left out and counted in the record.
- `rDNA.fa`: the NCBI rDNA unit (human U13369.1, mouse BK000964.3, rat NR_046239.1) and the Ensembl rRNA and Mt_rRNA genes. `tDNA.fa`: the GtRNAdb tRNA genes taken from the genome; a tRNA whose sequence differs from GtRNAdb's stops the build (wrong assembly).

The genome name is `gcf`, not one of TEsmall's built-in names, so a database that fails TEsmall's own checks stops the run instead of starting a download. TEsmall removes the string `genomes` from its database path, so an `EXT_DIR` containing it stops the workflow at load.

**Run** (`tesmall_quant`, one job per sample in its own directory, as TEsmall writes into its working directory and always names its table `count_summary.txt`): `gcfntnu/tesmall:2.1.0-1` with `--no-trim` ([gcfntnu/gcf-docker#47](https://github.com/gcfntnu/gcf-docker/pull/47)) on the shared filtered reads. Unpatched TEsmall trims every input again, including a hard-coded TruSeq 5′ adapter whose partial matches it would cut from the start of already-trimmed reads (49 human mature miRNAs begin with one); with `--no-trim` only its length filter remains: the kit's filtered minimum to `quant.tesmall.max_length`. TEsmall cuts sample labels at dots, so such sample names stop the workflow at load. It runs single-threaded.

**Outputs** in `quant/tesmall/`: `<sample>/count_summary.txt` with TEsmall's other per-sample files, and `tesmall_counts.tsv` (`ftype`, `fid`, one column per sample; [`tesmall_counts.py`](rules/quant/scripts/tesmall_counts.py)).

On the benchmark's reduced reference (mouse chromosome 19) the database build took under a minute, and placed all 24 piRNAclusterDB clusters of chromosome 19 (published on GRCm38) on GRCm39. On its synthetic reads TEsmall counted most miRNAs exactly; it gave the reads of `mmu-miR-3084-3p`, whose two loci lie inside snoRNA genes, to those snoRNAs (structural RNA comes first), and left a 3′ isomiR with a non-templated A unaligned (exact genome alignment). The full-genome database (a Bowtie index of the whole genome) has not been built yet.

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | `unitas`. Accepts `unitas`, `mirge3`, `mirflowz` and `tesmall`, alone or comma-separated; selects the `quant_all` targets only, not the delivered UNITAS outputs | [smallrna.config](smallrna.config) |
| `quant.mirge3.db` | `mirbase`. `mirbase` or `mirgenedb` (any case), passed to miRge3 as `miRBase` or `MirGeneDB`; other values stop the workflow at load | `smallrna.config` |
| `db.mirge3.release` | `09-11-2026` (miRBase 23.1, MirGeneDB 3.0, Ensembl 116); `sourceforge` for the old unversioned library. Unset stops the workflow at load | [main.config](../main.config) |
| `db.mirflowz.ensembl_release` | Per organism: `homo_sapiens` and `mus_musculus` `116`, `rattus_norvegicus` `104`; must be on the miRBase GFF's assembly, checked at load when mirflowz is selected | [main.config](../main.config) |
| `quant.mirflowz.mem_gb` | `60`. Memory budget (GB) of mirflowz's nested runs, preparation and quantification; per-rule memory and the total of concurrent jobs are capped at it | `smallrna.config` |
| `db.mirflowz.custom` | `name`, `genome`, `gtf`, `mirna_gff`, all `null`; set together to replace the mirflowz reference inputs (also for the TEsmall database) | [main.config](../main.config) |
| `quant.tesmall.max_length`, `mismatch`, `maxaln` | `36`, `0`, `100`: TEsmall's maximum read length (the minimum is the kit's filtered minimum), genome mismatches (0–3) and the most alignments a read may have | `smallrna.config` |
| `align.star.args`, `align.star.guided_index` | STAR options for short reads | `smallrna.config` |
| `align.genome.aligner` | `star`; any other value stops the workflow at load | [smallrna.config](smallrna.config) |
| `filter.dedup` | `false`; `true` stops the workflow at load (no read deduplication) | `main.config`, kit entry |
| `filter.trim` | Per kit | Kit entry ([#215](https://github.com/gcfntnu/gcf-workflows/issues/215)) |
| `filter.spikein.ref` | `null` (no spike-ins, no spike-in step). The one key a project with spike-ins sets. Otherwise comma-separated spike-in sets, each `ERCC92`, `smallrna_calibrators` or the path of a project FASTA (`.fa`, `.fasta`), e.g. `ERCC92,smallrna_calibrators`; the bowtie route, or UNITAS and miRge3 (`quantifier: native`), use these sets. Unknown sets, a set listed twice, an unreadable or invalid FASTA, or `quantifier: native` without any set stop the workflow at load | `smallrna.config` |
| `filter.spikein.quantifier` | Who estimates the spike-ins, used only with sets listed: `bowtie` (default), a prefilter that counts each spike-in and removes those reads before quantification ([spike-in counts](#spike-in-counts-bowtie)); `native`, every quantifier with spike-in support estimates them from the reads with the spike-ins (UNITAS with `-refseq`, miRge3 with a spike-in index), while the bowtie prefilter still counts them and removes them for every other step. `skip` (no spike-ins is now an empty `filter.spikein.ref`), `unitas` (the former name of `native`), `filter.spikein.filter` (now implied by `bowtie`) and other values stop the workflow at load | `smallrna.config` |
| `filter.trim.quantifier` | `cutadapt` (kit default, [filtering with cutadapt](#filtering-with-cutadapt-kit-default)), `fastp` ([fastp route](#fastp-route-filtertrimquantifier-fastp)) or `skip`; any other value stops the workflow at load (the miRTrace filter route was removed) | Kit entry |
| `filter.trim.cutadapt` | Per kit: `adapter5`, `five_prime_piece`, `random_ends`, `minimum_length`, `max_n`, `max_ee`; required with `filter.trim.quantifier: cutadapt` or `skip`, set for every small-RNA kit | Kit entry |
| `filter.trim.fastp.args` | No kit default; required with `filter.trim.quantifier: fastp`, and set with another route it stops the workflow at load | Project |
| `filter.contaminants` | `skip`/off; `quantifier: unitas` adds the sequences to the UNITAS reference | `smallrna.config` |
| `filter.ribosomal.quantifier` | `skip`, the only accepted value. UNITAS no longer supports `-riborase`, so `unitas` stops the workflow at load | `smallrna.config` |
| `db.reference_db`, `db.mirbase.release` | `ensembl`, releases in `main.config`. `db.mirbase.release` is `23.1`, the release the `05-11-2026` UNITAS dumps use, and must equal the selected miRBase dump's release (checked by `unitas_reference`); miRBase also serves `20`, `21` and `22`, and other values stop the workflow at load | Kit entry, `main.config` |
| `db.unitas.mirna_db` | `mirbase`. `mirgenedb` uses a UNITAS dump built from MirGeneDB instead of miRBase (see below); other values stop the workflow at load | [main.config](../main.config) |
| `db.unitas.release` | `05-11-2026`; with `mirna_db: mirbase`, the reference dump `UNITAS_refdump_<organism>` is downloaded from this release | [main.config](../main.config) |
| `db.unitas.mirgenedb_release` | `06-11-2026-mirgenedb` (MirGeneDB 3.0); with `mirna_db: mirgenedb`, the dump is downloaded from this release. It must hold MirGeneDB dumps (`gcf_mirna_db.tsv`), and unset stops the workflow at load | [main.config](../main.config) |

Precedence between these sources is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `smallrna/.tests`: synthetic dry-run fixture (NEXTflex v3, `mus_musculus`, four samples, with small reference FASTA stubs). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-small-rna.yml`).
- `smallrna/.tests/test_fastp.py`: the fastp route with the kits' former presets: configuration, all-kit dry runs, stage accounting, exact v3 clipping and 15/16/17 nt boundaries, quality/N-base/complexity and polyG outcomes with fastp 0.20.1. All inputs are synthetic; no new v4 fixture is built.
- `smallrna/.tests/test_cutadapt.py`: the cutadapt route's settings and their validation, commands per kit, the truncated 5′ adapter set, the merged report and the route's dry runs (the default for every kit; projects with fastp settings keep fastp; fastp settings with another route stop at load); with cutadapt 5.2 on synthetic reads, that 5′-adapter fragments and adapter dimers are removed while miR-22-5p and miR-139-5p are kept, and that a 3′ adapter with a deletion is trimmed.
- `smallrna/.tests/test_run_unitas.py` and `test_bfq_unitas_links.py`: per-sample UNITAS work directories with a stand-in `unitas.pl` (concurrent and similarly named samples, interrupted runs, unexpected results, non-zero exit, `ERROR!` with exit 0, missing annotation summary), and relative delivery links.
- `smallrna/.tests/test_route_guards.py`: dry runs showing that a non-STAR aligner, `filter.dedup: true`, an organism without a UNITAS reference dump or `filter.ribosomal.quantifier: unitas` stops at load, while the defaults plan delivery; and that `-refseq` reaches UNITAS only when spike-ins or contaminants are configured; miRTrace reads each route's filtered reads (fastp, cutadapt, skip) without an adapter (`-p illumina`, two-column config file), and the MultiQC template shows only its contamination and complexity sections.
- `smallrna/.tests/test_kit_metadata.py`: read orientation of the kits, the NEXTflex v4 `blockers` default, its report-header argument and the load errors.
- `smallrna/.tests/test_spikein_bowtie.py`: the bowtie spike-in route: the combined reference (DNA, unique names), complete count tables with ambiguous reads, the matrix; dry runs of counting, removal for every consumer, project FASTA sets (DNA copy for every route, a path with a space) and the load errors; the MultiQC General Statistics totals for both routes; and, with bowtie 1.x on `PATH`, that exact forward pieces count while reads with a mismatch, reverse-complement or longer reads do not, and that a read matching one spike-in twice and another once is ambiguous.
- `smallrna/.tests/test_quant_methods.py`: `quant_all` dry runs for each supported `quant.method`, miRge3 on the filtered reads per kit (no own trimming or deduplication, the kit's minimum length) and the merged reads with trimming skipped, and errors for removed or unknown values; the dump verification and reference check planned before the UNITAS AnnData, for miRBase and MirGeneDB dumps.
- `smallrna/.tests/test_mirflowz_route.py`: mirflowz dry runs: combined route with delivery unchanged, reuse of the project's Ensembl files for mouse, Ensembl 104 downloads and the removed contig for rat, custom references, kit read settings, and load errors only when mirflowz is selected.
- `smallrna/.tests/test_tesmall_route.py`: TEsmall dry runs: the database on the mirflowz reference with its GtRNAdb, RepeatMasker, rDNA and UNITAS piRNA cluster sources (mouse and rat), `--no-trim` and the run settings, delivery unchanged, and the load errors (organism, maximum length, dotted sample names); the joined count table.
- `gcfdb/.tests/test_tesmall_db.py`: the TEsmall database classes from synthetic miRBase, GTF, RepeatMasker and GtRNAdb inputs, piRNA cluster headers and placement from anchor hits, and, with samtools and Bowtie 1 on `PATH` (e.g. in the TEsmall image), a build on a synthetic genome including the tRNA sequence check.
- `smallrna/.tests/test_run_mirflowz.py` and `gcfdb/.tests/test_mirflowz_prepare.py`: the nested mirflowz runs with a stand-in `snakemake`: only the needed rule files included, linked prepared reference, read settings, input conversion, outputs kept, missing files and failures reported with the tool log.
- `smallrna/.tests/test_benchmark.py`: the quantification benchmark ([benchmark/](benchmark/README.md)): synthetic truth, method adapters, metrics, and prepare/evaluate on synthetic tables.
- `gcfdb/.tests/test_mirflowz_chr_map.py`: chromosome maps from synthetic genome and GFF names: UCSC to Ensembl names, unplaced contigs, unmapped names, explicit removal and assembly mismatch.
- `smallrna/.tests/test_mirge3_counts.py`: miRge3 counts from a synthetic `mapped.csv`: exact and isomiR reads together, miRNAs kept with 2 exact reads in the project (also 1 per sample) and dropped with 1 or only isomiRs, merged names and SNP variants folded, and errors for reads with two names, unknown miRNAs or differing samples.
- `smallrna/.tests/test_unitas_modifications.py`: the four UNITAS modification tables are written also when no sample has a modification type (header only), and tables with rows are unchanged.
- `smallrna/.tests/test_unitas_reference.py`: dump records against synthetic annotations: matching miRBase and MirGeneDB releases, release and database mismatches, inconsistent or unreadable dump records, and missing records reported as `NA`.
- `gcfdb/.tests/test_unitas_verify.py`: a cached dump against a synthetic published archive: verified, checksum mismatch, edited, extra or missing files, no `SHA256SUMS`, archive not listed or not downloadable.
- `smallrna/.tests/test_create_anndata.py`: AnnData regression, run in CI on Python 3.14 (see below).
- `smallrna/.tests/test_report_figures.py`: report figures: the PCA on log2(CPM + 1) (a library differing only in depth lands on its original; explained variance on the axes), the top-miRNA plot in counts per million miRNA reads, the components from the filtered matrix, a note when two components are impossible (two samples; no miRNA left), and the PCA section rewritten for every sample count, also over an earlier four-sample PCA (dry runs). The figure tests run in CI with the scanpy image's scanpy and anndata versions.
- Commands and environments: [development guide](../docs/development.md#run-the-checks).

## Report figures

`bfq_level2_pca` writes the report's PCA from the UNITAS miRNA counts: log2(CPM + 1), with CPM over all miRNA reads, of the miRNAs detected in at least min(3, samples) samples; the axes give the share of variance explained. The number of components comes from the filtered matrix; with fewer than three samples, or fewer than three miRNAs detected in at least three samples, the section is a short note instead of the plot, so the report still builds and an earlier PCA is replaced. `bfq_level2_mirna_high` shows the ten miRNAs with the highest mean counts per million miRNA reads, as CPM and as read counts. Neither uses the DESeq2 layer below.

## UNITAS AnnData and variance transformation

The UNITAS workflow writes raw counts as float32 in `quant/unitas/adata.h5ad`, with samples in `obs` and miRNA features in `var`. Sample and feature metadata are aligned to the count table; absent feature annotations remain missing and emit a warning rather than dropping count rows. H5AD metadata column names escape `%` as `%25` and `/` as `%2F` (for example `260/280` becomes `260%2F280`). Original names are retained in `uns['metadata_column_names']['obs' or 'var']`, keyed by the encoded name. Cell values and sample/feature identifiers are unchanged.

`unitas_mirna_vsn` computes DESeq2-transformed counts in the existing `tximport` image and writes `mirna_vsd.tsv`. The AnnData job loads this table into `layers['vsd']` in the Scanpy image, which does not need R or rpy2. Both axes must exactly match the count table. The transformation uses features by samples, `round(counts + 1)` and `fitType="mean"`; it chooses rlog for at most 50 samples and variance-stabilizing transformation for larger cohorts. `uns['varians_transform']` records the method. This corrects the old in-process conversion's reversed feature/sample axes without changing raw counts.

The standalone creator still supports `--add-vsn` with R, DESeq2 and rpy2 installed. Use `--add-vsn --vsn-input transformed.tsv` for a precomputed table produced by `mirna_vsn.R` with the same settings. The existing versioned Scanpy and tximport image selections are retained; the API regression job tests AnnData 0.13.3.post0 with NumPy 2.5.2 and pandas 3.0.5 rather than downgrading the shared image.

Run the Python regressions with:

```bash
python -m pytest -q smallrna/.tests/test_create_anndata.py
```

These tests cover CLI/H5AD roundtrips, count dtype, metadata ordering and VSN layer assembly. The standalone R bridge is mocked in these tests; real DESeq2 and complete human/mouse workflows require container integration.

### Credits

## More documentation

- [`methods.txt`](methods.txt): methods text for reports. It is not read by the workflow.
- [Architecture](../docs/architecture.md), [configuration](../docs/configuration.md), [reference data](../gcfdb/README.md)

## Credits

Written by Arnar Flatberg ([flatberg](https://github.com/flatberg)) at the Genomics Core Facility, [Norwegian University of Science and Technology, NTNU](http://www.ntnu.no). Licence: [LICENSE](LICENSE).
