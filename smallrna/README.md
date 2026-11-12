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

The kit entries set `filter.trim.fastp.args`, miRTrace `protocol: nextflex` and `db.reference_db: ensembl`.

### fastp stages and final insert policy

All three kits use **adapter removal, then final trimming/filtering**, with fastp 0.20.1.
The workflow generates the adapter-only options `-Q -L -G`: no quality/N-base filtering, optional length filtering or polyG trimming in that pass.
Empty/null reads are still rejected by fastp and counted in that stage's `too_short_reads`; `-L` does not preserve empty reads.
Only the second invocation receives `filter.trim.fastp.args`; the workflow supplies `-A` there to prevent repeated adapter removal.
There is no configurable first-stage option string, and no small-RNA use of the shared rule's `params` key.

| Kit | Fixed clipping after adapter removal | Minimum final insert |
| --- | --- | --- |
| NEXTflex v3 SE | Four bases from each exposed end (`-f 4 -t 4`) | 16 nt |
| NEXTflex v4 SE | None | 15 nt (existing policy preserved) |
| TruSeq Small RNA SE | None | 15 nt (existing policy preserved) |

Every preset explicitly sets the existing final quality/N-base and low-complexity thresholds:
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

The [v4 manufacturer protocol](https://resources.revvity.com/pdfs/man-nextflex-small-rna-trimming-instructions.pdf) recommends cutadapt quality cutting at Q20 and a 16 nt minimum.
That is a possible future v4 policy, not the current preset, and Q20 end cutting is not equivalent to fastp's read-level quality filtering.
This change preserves v4's existing policy and does not claim protocol equivalence or require new v4 fixture data.
The v3 16 nt minimum is the workflow requirement in [#215](https://github.com/gcfntnu/gcf-workflows/issues/215), not a claim about that v4 protocol.
References, alignment and quantification are outside this change.

### Existing project overrides

New project configurations receive the kit's final options under the existing string key `filter.trim.fastp.args`.
Projects retaining the old sibling `trim_args` (even empty), or a nonempty `params`, fail early with migration guidance instead of silently changing the stages.
For an intentional migration, inspect the old effective configuration and remove `trim_args`; do **not** append its `-Q -n 0 -l 0 -G` to the final options.
Those were first-pass controls. The workflow now supplies equivalent adapter-only controls independently of final filtering.
Set `args` to the applicable preset above, including v3's `-f 4 -t 4 -l 16`, or preserve explicitly reviewed project-specific final overrides.
The argument string must specify a final length policy (`-l N` or `-L`, not both); use `-l 0` or `-L` only for deliberately disabled minimum-length filtering.
In fastp, `-Q` also disables N-base filtering, regardless of `-n`; document that choice when overriding the final policy.
Resuming an older workdir uses its captured workflow/configuration; deploying this revision does not migrate it automatically.

### Reports

The existing `<sample>.fastp.json` delivery path now contains reconciled original-input/final-output metrics, both commands in execution order,
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
| filter | `rules/filter/cliptrim.smk`, shared `common/` rules | FASTQ merging; clipping and trimming |
| align | `rules/align/star.smk` | STAR (smallrna options); the only supported genome aligner |
| quant | `rules/quant/spikein.smk`, `unitas.smk`, `mirge3.smk` | Spike-in sets shared by both quantifiers; UNITAS (including its variance-transformation rules); miRge3 |
| qc | `rules/qc/`, shared `common/` rules | FastQC (shared); miRTrace QC on the merged FASTQs, delivered to `bfq/logs/mirtrace/` and shown in the MultiQC report |
| bfq | `rules/bfq/` (levels 1–3) | Collect delivery products into `bfq/` |

**`quant.method` and delivery.** `quant.method` selects what the separate `quant_all` target builds: `unitas`, `mirge3` or both (comma-separated). Any other value stops `quant_all` with an error. It does **not** change delivery: `bfq_all` and `multiqc_report` always deliver the UNITAS products (`rules/bfq/bfq_level2.smk`), whatever `quant.method` says.

**Inactive source.** Files present in the source tree are not necessarily loaded. The entry point does not load `rules/filter/calibrators.smk`.

**Quantifier status.** UNITAS is no longer distributed upstream, and its download server is gone. The workflow uses `gcfntnu/unitas:1.8.0-1`, UNITAS 1.8.0 patched to run offline and to exit non-zero on errors, with reference dumps built from Ensembl, GtRNAdb and miRBase (`db.unitas.release`). Dumps exist only for `homo_sapiens`, `mus_musculus` and `rattus_norvegicus`; any other organism stops the workflow at load. UNITAS remains the delivery route until a maintained quantifier replaces it. miRge3 0.1.4 is kept as an alternative `quant_all` route, with versioned libraries (`db.mirge3.release`, see [miRge3 libraries](#mirge3-libraries)). miRge (v2), QuickMIRSeq and SeqBuster were removed.

**UNITAS execution.** Each sample runs in its own work directory, `quant/unitas/work/<sample>/`, with the reference dump linked in ([`run_unitas.sh`](rules/quant/scripts/run_unitas.sh)); `seqmap.exe` comes from the image. Exactly one UNITAS result directory must appear there; it is moved to `quant/unitas/unitas/<sample>/` and the work directory is removed. Concurrent or similarly named samples and leftovers from interrupted runs therefore cannot mix. UNITAS console output goes to `quant/unitas/logs/<sample>.log`. The job fails, and prints the end of that log, when UNITAS exits non-zero, reports `ERROR!` or writes no `unitas.annotation_summary.txt`. Spike-in and contaminant sequences are passed to UNITAS with `-refseq` only when `filter.spikein.quantifier` or `filter.contaminants.quantifier` is `unitas`.

## Outputs in `bfq/`

`data/tmp/smallrna/bfq/`, as planned by the fixture's dry run (NEXTflex v3). The set is the same for every `quant.method`, because delivery always uses UNITAS.

```text
multiqc_<project>.html, .multiqc_config.yaml
align/<sample>.sorted.bam
exprs/mir_counts.tsv, isomir_counts.tsv, isomir_counts_anno.tsv, trf_counts.tsv, isotrf_counts.tsv, annotations.tsv
figs/gene_high_mqc.yaml, pca_mqc.yaml
logs/unitas/<sample>  (relative links to quant/unitas/unitas/<sample>)
provenance/fastp/<sample>.adapter.{json,html}, <sample>.filter.{json,html}
logs/mirtrace/mirtrace-results.json, mirtrace-report.html, mirtrace-stats-{length,contamination_basic,mirna-complexity}.tsv
logs/<sample>/<sample>.fastp.html, .fastp.json, _fastqc.html, _fastqc.zip, .Log.final.out
```

**MirGeneDB dumps.** With `db.unitas.mirna_db: mirgenedb`, UNITAS annotates miRNAs against MirGeneDB 3.0 instead of miRBase: the curated set (human 514 precursors, 559 mature and 468 star sequences, against about 1,900 human miRBase hairpins), with the other 113 MirGeneDB species for cross-species hits. The rest of the dump (Ensembl, GtRNAdb, piRNA clusters, tRF-1, tRNA leaders, SILVA) is the same as in the miRBase dump. miRNAs are named by MirGeneDB gene and arm, without the species code: `Let-7-P1b-5p` and `Let-7-P1b-3p` for MirGeneDB's `Hsa-Let-7-P1b_5p` and `Hsa-Let-7-P1b_3p*`. That is the name UNITAS also gives reads it places on a hairpin arm (isomiRs), so each arm is counted under one name. Tables are not comparable by name with miRBase runs. `mirgenedb_genes` annotates the AnnData from the dump's MirGeneDB GFF: `star` marks the star arm, and `alias` holds the miRBase accession where there is one. UNITAS' `db_versions.info` still calls the database miRBase; the UNITAS log records `UNITAS MirGeneDB,3.0`. miRge3 keeps using `quant.mirge3.db`.

**UNITAS reference checks.** Two rules identify the selected dump and keep its counts joined to a matching annotation; neither downloads or rebuilds a dump already in `EXT_DIR`, and `unitas_db` itself is unchanged.

- `unitas_db_verify` (once per dump, beside it in `EXT_DIR/unitas/<release>/verify/UNITAS_refdump_<organism>.tsv`): when the release publishes `SHA256SUMS` (releases from `04-11-2026` on), the archive is downloaded, its SHA-256 must equal the published one, and the cached dump must hold exactly the archive's files with the same contents. A difference stops the workflow and names the files; the cache is left as it is. Releases without `SHA256SUMS` are recorded as `not verified`, and remain usable.
- `unitas_reference` (per project, `quant/unitas/reference.tsv`): records the dump's miRNA database and release (from `gcf_mirna_db.tsv` where the dump has it, otherwise UNITAS' `db_versions.info`), its Ensembl, GtRNAdb, piRNA cluster, tRF-1, tRNA leader and SILVA versions, whether it has source checksums (`gcf_sources.tsv`), the verification status, and the annotation's release, each with the file it comes from. Records a dump does not have are `NA`; nothing is filled in from current settings. It is an input of `unitas_mirna_anndata` and `unitas_all`.

UNITAS miRNA counts are joined to the feature annotation by name, and names change between releases. On synthetic reads of every human mature miRNA, the `05-11-2026` dump (miRBase 23.1) leaves 6 of 2,730 counted miRNAs unannotated with the miRBase 23.1 GFF (mature sequences without a genome locus there), but 963 (35%) with the miRBase 22 GFF; MirGeneDB counts match none of the miRBase names. `unitas_reference` therefore stops when the dump's miRNA database is not `db.unitas.mirna_db`, when a miRBase dump's release differs from `db.mirbase.release` (the error names both, and suggests `db.mirbase.release` only for releases gcfdb can download), when a MirGeneDB dump disagrees with its own GFF, or when the dump's miRNA release cannot be read. The dump's Ensembl release is recorded, not checked: UNITAS' reference need not match the genome alignment reference.

## miRge3 libraries

`db.mirge3.release` selects the miRge3 library; `quant.mirge3.db` (`mirbase` or `mirgenedb`) selects the miRNA database within it. The default release, `09-11-2026`, holds libraries for human, mouse and rat built by [gcfntnu/gcf-docker#41](https://github.com/gcfntnu/gcf-docker/pull/41) from miRBase 23.1, MirGeneDB 3.0 and Ensembl 116, in the format of the SourceForge libraries. Each library records its sources: `gcf_mirna_db.tsv` (miRNA releases), `gcf_manifest.tsv` (releases, assemblies, tools) and `gcf_sources.tsv` (every download with URL and SHA-256).

- `mirge3_library` (once per release and organism, [`gcfdb/scripts/mirge3_library.py`](../gcfdb/scripts/mirge3_library.py)): downloads `miRge3_lib_<organism>.tar.gz` from `<winecellar.url>/mirge3/<release>/`. The archive's SHA-256 must equal the release's `SHA256SUMS` before it is unpacked, and the library must hold every file miRge3 opens for both databases (and the tRNA fragment annotation for human), with all six files of every Bowtie index, plus the three provenance files. It is then installed as `EXT_DIR/mirge3/<release>/<org>/` and recorded in `EXT_DIR/mirge3/<release>/verify/miRge3_lib_<organism>.tsv`. A missing `SHA256SUMS`, a checksum or file mismatch, or an existing `<org>/` without that record stops the workflow; each release has its own directory, so other caches are not touched.
- `db.mirge3.release: sourceforge` selects the unversioned SourceForge library at its old place, `EXT_DIR/mirge3/<organism>/<org>/`, for retained analyses. Existing copies are used as they are; the library has no checksums or provenance.
- `mirge3_reference` (per project, `quant/mirge3/reference.tsv`): stops when the library lacks the files of `quant.mirge3.db` (including any of the six files of its Bowtie indexes), or when its `gcf_mirna_db.tsv` does not list that database. It records the database and its release, the Ensembl release and assemblies from the manifest, the source checksums and the archive check, each with the file it comes from; the SourceForge library's are `NA`. miRge3 counts are not joined to another annotation, so no other release is checked. It is an input of `mirge3_quant` and `quant_all`.

**Spike-ins.** Spike-ins are a library prep choice: `filter.spikein.ref` lists the sets, and their sequences come from gcfdb (`EXT_DIR/spikein/fasta/<set>.fa`) for either quantifier. miRge3 0.1.4 counts spike-ins (`--spikeIn`) only from the library's `index.Libs/<org>_spike-in` index and takes no sequences itself. The published library's spike-in index is therefore removed when it is installed, and the shared library is not changed after that. Each project gets its own library, `quant/mirge3/lib/<org>/`, made by `mirge3_project_library` ([`mirge3_project_library.py`](rules/quant/scripts/mirge3_project_library.py)): links to every file of the shared library, plus a spike-in index built from the project's sets, or none when no sets are listed; `--spikeIn` is passed only then, and miRge3's `--libraries-path` is `quant/mirge3/lib`. Projects therefore cannot use each other's spike-ins, and a change of sets rebuilds the project library. `quant/mirge3/spikein.tsv` records the shared library, the sets, their sequence count and the SHA-256 of the combined FASTA. miRge3 aligns spike-ins last, exactly, and only for reads nothing else matched, and reports one `Spike-in` total per sample in `annotation.report.csv`; counts per spike-in come from UNITAS (`filter.spikein.quantifier: unitas`).

**Changes from the SourceForge libraries.** The new libraries are a scientific reference update; counts are not comparable by name with SourceForge runs, and unchanged counts are not expected. From gcf-docker#41's comparison:

- miRBase 22 to 23.1: renames (`hsa-miR-1-3p` to `hsa-miR-1-2-3p`, `hsa-miR-107` to `hsa-miR-107-3p`) and mature boundaries moved by 1-2 nt (`mmu-let-7g-5p` gains a 3' T). Same name with changed sequence / renamed / dropped / added: human 451 / 743 / 120 / 252, mouse 373 / 269 / 92 / 93, rat 220 / 117 / 30 / 132.
- Merge groups are miRNAs with identical mature sequences and a miRNA's SNP variants (for `--AtoI`, from dbSNP). The SourceForge libraries also merged similar but not identical sequences by an undocumented rule; such miRNAs (e.g. `miR-320a/b/c/d/e`) are now reported separately. Human miRBase merge groups: 284 to 226.
- MirGeneDB 2.0 to 3.0: renamed families (`Hsa-Let-7-P2_5p` to `Hsa-Let-7-P1b_5p`), and the `-v1`/`-v2` arm variants are gone (human mature 1171 to 1027).
- Other RNA from Ensembl 116 (human mRNA 170,448 to 465,253 sequences; `rRNA_pseudogene` counted as rRNA, as in the old library). Mouse and rat tRNAs from GtRNAdb; the human tRNA set and its tRNA fragment files are unchanged. The old rat library was far older (Rnor_5.0 transcripts).
- No genome pickles, so `--novel-miRNA` is not available (not used by the workflow).
- miRge3 0.1.4 writes `miRBase22` or `MirGeneDB2.0` into its GFF output whatever the library; `reference.tsv` records the actual release.

**Delivered UNITAS tables.** `mir_counts.tsv`, `trf_counts.tsv` and `isomir_counts.tsv` have features as rows and samples as columns; `isomir_counts.tsv` rows are `sequence|miR-name`, annotated in `isomir_counts_anno.tsv`. `isotrf_counts.tsv` (tab-separated) and `annotations.tsv` are long tables with a `Sample_ID` column. Every `annotations.tsv` row names its parent levels (`Level1`, `Level2`, `Level3`). All tables are sorted, so repeated runs give identical files.

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | `unitas`. Accepts `unitas`, `mirge3` or both, comma-separated; selects the `quant_all` targets only, not the delivered UNITAS outputs | [smallrna.config](smallrna.config) |
| `quant.mirge3.db` | `mirbase`. `mirbase` or `mirgenedb` (any case), passed to miRge3 as `miRBase` or `MirGeneDB`; other values stop the workflow at load | `smallrna.config` |
| `db.mirge3.release` | `09-11-2026` (miRBase 23.1, MirGeneDB 3.0, Ensembl 116); `sourceforge` for the old unversioned library. Unset stops the workflow at load | [main.config](../main.config) |
| `quant.mirge3.uniq_mol_ids` | `null`; NEXTflex v3 sets `'4,4'` (4N randomized ends). Passed to miRge3 as `--uniq-mol-ids` | `smallrna.config`, kit entry |
| `align.star.args`, `align.star.guided_index` | STAR options for short reads | `smallrna.config` |
| `align.genome.aligner` | `star`; any other value stops the workflow at load | [smallrna.config](smallrna.config) |
| `filter.dedup` | `false`; `true` stops the workflow at load (no read deduplication) | `main.config`, kit entry |
| `filter.trim` | Per kit | Kit entry ([#215](https://github.com/gcfntnu/gcf-workflows/issues/215)) |
| `filter.spikein.ref` | `null` (no spike-ins). Otherwise comma-separated spike-in sets, each `ERCC92` or `smallrna_calibrators`, e.g. `ERCC92,smallrna_calibrators`; UNITAS (with `quantifier: unitas`) and miRge3 use these sets. Unknown sets, or `filter.spikein.quantifier: unitas` without any set, stop the workflow at load | `smallrna.config` |
| `filter.trim.quantifier` | `fastp` or `skip`; any other value stops the workflow at load (the miRTrace filter route was removed) | Kit entry |
| `filter.trim.mirtrace.protocol` | Per kit: `nextflex` (NEXTflex v3, 4N ends), `illumina` (NEXTflex v4, TruSeq). Required | Kit entry |
| `filter.spikein`, `filter.contaminants` | `skip`/off; `quantifier: unitas` adds the sequences to the UNITAS reference | `smallrna.config` |
| `filter.ribosomal.quantifier` | `skip`, the only accepted value. UNITAS no longer supports `-riborase`, so `unitas` stops the workflow at load | `smallrna.config` |
| `db.reference_db`, `db.mirbase.release` | `ensembl`, releases in `main.config`. `db.mirbase.release` is `23.1`, the release the `05-11-2026` UNITAS dumps use, and must equal the selected miRBase dump's release (checked by `unitas_reference`); miRBase also serves `20`, `21` and `22`, and other values stop the workflow at load | Kit entry, `main.config` |
| `db.unitas.mirna_db` | `mirbase`. `mirgenedb` uses a UNITAS dump built from MirGeneDB instead of miRBase (see below); other values stop the workflow at load | [main.config](../main.config) |
| `db.unitas.release` | `05-11-2026`; with `mirna_db: mirbase`, the reference dump `UNITAS_refdump_<organism>` is downloaded from this release | [main.config](../main.config) |
| `db.unitas.mirgenedb_release` | `06-11-2026-mirgenedb` (MirGeneDB 3.0); with `mirna_db: mirgenedb`, the dump is downloaded from this release. It must hold MirGeneDB dumps (`gcf_mirna_db.tsv`), and unset stops the workflow at load | [main.config](../main.config) |

Precedence between these sources is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `smallrna/.tests`: synthetic dry-run fixture (NEXTflex v3, `mus_musculus`, four samples, with small reference FASTA stubs). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-small-rna.yml`).
- `smallrna/.tests/test_fastp.py`: configuration, all-kit dry runs, stage accounting, exact v3 clipping and 15/16/17 nt boundaries, quality/N-base/complexity and polyG outcomes with fastp 0.20.1. All inputs are synthetic; no new v4 fixture is built.
- `smallrna/.tests/test_run_unitas.py` and `test_bfq_unitas_links.py`: per-sample UNITAS work directories with a stand-in `unitas.pl` (concurrent and similarly named samples, interrupted runs, unexpected results, non-zero exit, `ERROR!` with exit 0, missing annotation summary), and relative delivery links.
- `smallrna/.tests/test_route_guards.py`: dry runs showing that a non-STAR aligner, `filter.dedup: true`, an organism without a UNITAS reference dump or `filter.ribosomal.quantifier: unitas` stops at load, while the defaults plan delivery; and that `-refseq` reaches UNITAS only when spike-ins or contaminants are configured.
- `smallrna/.tests/test_quant_methods.py`: `quant_all` dry runs for each supported `quant.method`, and errors for removed or unknown values; the dump verification and reference check planned before the UNITAS AnnData, for miRBase and MirGeneDB dumps.
- `smallrna/.tests/test_unitas_reference.py`: dump records against synthetic annotations: matching miRBase and MirGeneDB releases, release and database mismatches, inconsistent or unreadable dump records, and missing records reported as `NA`.
- `gcfdb/.tests/test_unitas_verify.py`: a cached dump against a synthetic published archive: verified, checksum mismatch, edited, extra or missing files, no `SHA256SUMS`, archive not listed or not downloadable.
- `smallrna/.tests/test_create_anndata.py`: AnnData regression, run in CI on Python 3.14 (see below).
- Commands and environments: [development guide](../docs/development.md#run-the-checks).

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
