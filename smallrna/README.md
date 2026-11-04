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
| quant | `rules/quant/unitas.smk`, `mirge3.smk` | UNITAS (including its spike-in and variance-transformation rules), miRge3 |
| qc | `rules/qc/`, shared `common/` rules | FastQC (shared); miRTrace QC on the merged FASTQs, delivered to `bfq/logs/mirtrace/` and shown in the MultiQC report |
| bfq | `rules/bfq/` (levels 1–3) | Collect delivery products into `bfq/` |

**`quant.method` and delivery.** `quant.method` selects what the separate `quant_all` target builds: `unitas`, `mirge3` or both (comma-separated). Any other value stops `quant_all` with an error. It does **not** change delivery: `bfq_all` and `multiqc_report` always deliver the UNITAS products (`rules/bfq/bfq_level2.smk`), whatever `quant.method` says.

**Inactive source.** Files present in the source tree are not necessarily loaded. The entry point does not load `rules/filter/calibrators.smk`.

**Quantifier status.** UNITAS is no longer distributed upstream; the workflow uses an archived image and reference dump, so its annotations cannot be updated. It remains the delivery route until a maintained quantifier replaces it. miRge3 is kept as an alternative `quant_all` route. miRge (v2), QuickMIRSeq and SeqBuster were removed.

**UNITAS execution.** Each sample runs in its own work directory, `quant/unitas/work/<sample>/`, with `seqmap.exe` and the reference dump linked in ([`run_unitas.sh`](rules/quant/scripts/run_unitas.sh)). Exactly one UNITAS result directory must appear there; it is moved to `quant/unitas/unitas/<sample>/` and the work directory is removed. Concurrent or similarly named samples and leftovers from interrupted runs therefore cannot mix.

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

**Delivered UNITAS tables.** `mir_counts.tsv`, `trf_counts.tsv` and `isomir_counts.tsv` have features as rows and samples as columns; `isomir_counts.tsv` rows are `sequence|miR-name`, annotated in `isomir_counts_anno.tsv`. `isotrf_counts.tsv` (tab-separated) and `annotations.tsv` are long tables with a `Sample_ID` column. Every `annotations.tsv` row names its parent levels (`Level1`, `Level2`, `Level3`). All tables are sorted, so repeated runs give identical files.

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | `unitas`. Accepts `unitas`, `mirge3` or both, comma-separated; selects the `quant_all` targets only, not the delivered UNITAS outputs | [smallrna.config](smallrna.config) |
| `quant.mirge3.db` | `mirbase` | `smallrna.config` |
| `quant.mirge3.uniq_mol_ids` | `null`; NEXTflex v3 sets `'4,4'` (4N randomized ends). Passed to miRge3 as `--uniq-mol-ids` | `smallrna.config`, kit entry |
| `align.star.args`, `align.star.guided_index` | STAR options for short reads | `smallrna.config` |
| `align.genome.aligner` | `star`; any other value stops the workflow at load | [smallrna.config](smallrna.config) |
| `filter.dedup` | `false`; `true` stops the workflow at load (no read deduplication) | `main.config`, kit entry |
| `filter.trim` | Per kit | Kit entry ([#215](https://github.com/gcfntnu/gcf-workflows/issues/215)) |
| `filter.spikein.ref` | `null` (no spike-ins). Otherwise comma-separated spike-in sets, each `ERCC92` or `smallrna_calibrators`, e.g. `ERCC92,smallrna_calibrators`. Unknown sets, or `filter.spikein.quantifier: unitas` without any set, stop the workflow at load | `smallrna.config` |
| `filter.trim.quantifier` | `fastp` or `skip`; any other value stops the workflow at load (the miRTrace filter route was removed) | Kit entry |
| `filter.trim.mirtrace.protocol` | Per kit: `nextflex` (NEXTflex v3, 4N ends), `illumina` (NEXTflex v4, TruSeq). Required | Kit entry |
| `filter.spikein`, `filter.contaminants`, `filter.ribosomal` | `skip`/off | `smallrna.config` |
| `db.reference_db`, `db.mirbase.release` | `ensembl`, releases in `main.config` | Kit entry, `main.config` |

Precedence between these sources is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `smallrna/.tests`: synthetic dry-run fixture (NEXTflex v3, `mus_musculus`, four samples, with small reference FASTA stubs). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-small-rna.yml`).
- `smallrna/.tests/test_fastp.py`: configuration, all-kit dry runs, stage accounting, exact v3 clipping and 15/16/17 nt boundaries, quality/N-base/complexity and polyG outcomes with fastp 0.20.1. All inputs are synthetic; no new v4 fixture is built.
- `smallrna/.tests/test_run_unitas.py` and `test_bfq_unitas_links.py`: per-sample UNITAS work directories with a stand-in `unitas.pl` (concurrent and similarly named samples, interrupted runs, unexpected results), and relative delivery links.
- `smallrna/.tests/test_route_guards.py`: dry runs showing that a non-STAR aligner or `filter.dedup: true` stops at load, while the defaults plan delivery.
- `smallrna/.tests/test_quant_methods.py`: `quant_all` dry runs for each supported `quant.method`, and errors for removed or unknown values.
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
