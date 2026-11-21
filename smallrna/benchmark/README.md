# Small-RNA quantification benchmark

Compares the smallrna quantification methods (`quant.method`: UNITAS, miRge3, mirflowz, TEsmall, tiny-count,
nf-core/smrnaseq, biocore_htseq, biocore_miraligner) on datasets with a known truth, and makes adding a method or a
dataset a small change ([#338](https://github.com/gcfntnu/gcf-workflows/issues/338)).
It runs the workflow's own `quant_all` routes, never the tools directly, so it measures what a project would get.
It is not run in CI; its own code is tested in `smallrna/.tests/test_benchmark.py`.

## Steps

1. **Reference** (synthetic datasets). A reduced reference (a few chromosomes) lets mirflowz run on a workstation. It goes to
   `EXT_DIR/benchmark/references/<name>/`, the place the dataset declaration names:

   ```bash
   python smallrna/benchmark/reduce_reference.py \
       --genome https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.19.fa.gz \
       --gtf https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz \
       --mirna-gff https://www.mirbase.org/download_version_genome_files/23.1/mmu.gff3 \
       --mature https://www.mirbase.org/download_version_files/23.1/mature.fa \
       --species mmu --chromosomes 19 --output "$EXT_DIR/benchmark/references/mmu-chr19"
   ```

   UNITAS and miRge3 match reads by sequence and use their normal references from `EXT_DIR`.
2. **Project.** `benchmark.py prepare` writes a smallrna project from a declaration: configuration (organism, kit,
   `quant.method`, the reduced reference as `db.mirflowz.custom`), sample table, reads and the truth table:

   ```bash
   python smallrna/benchmark/benchmark.py prepare --dataset smallrna/benchmark/datasets/synthetic_mmu_chr19.yaml \
       --project <project> --ext-dir "$EXT_DIR"
   ```
3. **Run** `snakemake quant_all` in the project, as for any project on that machine (containers, cores, resources);
   with `--notemp` for pool and biological datasets, whose evaluation reads the workflow's temporary filtered reads
   (and, for biological datasets, the cutadapt reports; `evaluate` checks those before writing anything). After a
   complete run without `--notemp`, rerun with `--notemp --forceall`.
4. **Evaluate.** `benchmark.py evaluate --project <project>` writes `<project>/benchmark/`: `common.tsv` (every
   method's counts in one table), `isomirs.tsv` (the methods with isomiR output: their counts per read sequence, #454),
   `per_mirna.tsv` (pools: `per_sequence.tsv` and `filter_stages.tsv`), `metrics.tsv`
   and `report.md`; with UNITAS, `unitas_names.tsv` (how each UNITAS name was placed, below). UNITAS names are placed
   with the reference dump the project's UNITAS run recorded (`quant/unitas/reference.tsv`; read only; `--unitas-ref`
   names another).
   Before scoring, it stops when a method's count table does not have exactly the project's samples (a missing sample
   would count as zero), or, for pools and biological datasets, when the method's `quant/<method>/reference.tsv` names
   another miRNA database or release than the dataset's `mirbase.release` (#407). UNITAS and miRge3 record their
   reference's release; `quant_reference` records it for mirflowz, TEsmall, tiny-count and biocore_htseq (`db.mirbase.release`,
   or `NA` for a `db.mirflowz.custom` reference) and nf-core/smrnaseq and biocore_miraligner (`db.mirbase.release`), with the
   SHA-256 of the count
   table it was written after, which must match the table.
   Every metric gets a 95% bootstrap interval (#413), and so does every difference between two methods, paired
   because both are computed on the same resample. A resample draws the scoring units (truth miRNAs; the scored exact
   and shifted pool sequences) and, where there are several, the runs or samples, the ratiometric pools' runs within
   pool A and pool B, with replacement. `metrics_ci.tsv` (metric, row, value, low, high) and `differences.tsv`
   (metric, row_a, row_b, difference, low, high) are published with the metrics and listed with their SHA-256 in
   `evaluation.json`, with the resamples and seed (`--bootstrap N`, default 1000, `0` for none; `--seed`, default 1).
   `report.md` has a section with the intervals and the differences whose interval excludes 0. Some values have no
   interval: the synthetic totals `misassigned` and `unplaced`, and biological `detected_all_samples`, which a
   resample repeating a sample biases upwards. With few samples (e.g. 6 in GSE317995) the sample part of the interval
   is approximate.
   The last section of `report.md` is each method's runtime and memory (#386, `resources.py`), from the Snakemake
   benchmark files its rules write under `quant/benchmarks/`. Per method:
   - the wall time and CPU time of its quantification jobs, summed;
   - the largest job's peak memory;
   - the run's input reads per second of summed wall time.

   Reference steps built in the project, such as smrnaseq's hairpin index, are listed apart. Shared builds under
   `EXT_DIR` (the mirflowz reference, TEsmall's genome index and database, biocore_htseq's Bowtie 2 index) are not in the
   report: each records its build in its own folder, as `<reference>/benchmarks/<rule>.tsv`.
   Every job is in `resources.tsv`, which is listed in `evaluation.json`. A method run before #386, or with a step's
   file missing, is reported as not measured. `evaluate` checks the benchmark files before it writes anything.
   Snakemake rewrites a job's benchmark file whenever the job runs, so each file describes the run that made that
   job's current output. After a partial rerun a method's jobs can come from different runs, possibly on different
   machines; quantify a benchmark project from scratch to measure it in one run. The numbers depend on the machine, its storage and the cores and memory
   the run was given.
   The wall time is exact, I/O waits included. Snakemake samples memory and CPU every 30 s after a job's first 15 s,
   so short memory peaks can be missed. CPU time is a lower bound: a route that starts many short processes is
   undercounted most (mirflowz: 12 s recorded, 45 s measured with `/usr/bin/time`).

   **Differential expression** (#414). For a ratiometric pool that declares `de:`, `evaluate` also writes the inputs
   to `benchmark/de/input/`. With `de: {within: lab}` (#444) it also writes each lab's inputs to
   `benchmark/de/lab/<lab>/input/`, and lists the labs in `benchmark/de/levels.tsv`. `de.smk` tests every set it lists
   with edgeR's exact and GLM tests and limma-voom (`de.R`, in the `sartools` image), and `evaluate --de` adds the
   scores of their calls (Metrics, below):

   ```bash
   snakemake -s smallrna/benchmark/de.smk -d <project> --use-apptainer --cores 1
   python smallrna/benchmark/benchmark.py evaluate --project <project> --de
   ```

   Run `de.smk` with the same container settings as the project's `quant_all`: the repository holding `de.R` must be
   visible in the container. `de.R` fits copies of the inputs and records their MD5s, its own and that of the results
   table it writes. `evaluate --de` stops before writing anything unless that record matches the current inputs, `de.R`
   and the results table; after new counts, run `evaluate` without `--de` (it writes the current inputs), then
   `de.smk`, which also reruns after a change to `de.R`. A project prepared before #414 or #444 enables DE, or the
   tests within each lab, by adding its dataset's `de:` to `benchmark/dataset.yaml`. `de.smk` reads
   `benchmark/de/levels.tsv`, which `evaluate` writes since #444, so a project evaluated earlier needs `evaluate` once
   more before `de.smk`. The DE files are not listed in `evaluation.json`; the report gives each set's results SHA-256.
5. **Summary.** `benchmark.py summary` reads several evaluated projects and writes one overview, a Markdown file or a
   single self-contained HTML file (no external resources), by the output's extension:

   ```bash
   python smallrna/benchmark/benchmark.py summary --projects <project> [<project> ...] --output summary.html
   ```

   It has a section per dataset kind (synthetic, equimolar pool, ratiometric pool, biological) with a row per project:
   its dataset, kit, samples, reads processed, methods, a link to its
   `report.md`, and the kind's key metrics from its `metrics.tsv` per method (and the oracle, where it has values),
   each with its bootstrap interval in parentheses when the evaluation has one.
   A method not run on a dataset, or a metric missing from its table, is an empty cell. It writes nothing into the
   projects; an output inside one stops it. New evaluations publish `evaluation.json` with checksums for the metrics
   and report. For an older evaluation, the metrics must match the complete metrics table embedded in the report.
   New evaluations coordinate readers and writers through a regular lock file in the benchmark directory. Legacy
   summaries can read the existing report and metrics from a read-only project.
   The command prints the canonical summary path to open, including when its directory is a symlink.
   Rerun `evaluate` on the existing outputs if the pair is incomplete or inconsistent. `reads processed` is the reads
   the evaluated run filtered: each sample's input to the first cutadapt pass, from the JSON report the run wrote
   (`filter/cutadapt_stages/<sample>.adapter.json`, kept with `--notemp`), summed over the samples. `evaluate` reads
   them and records the total in `evaluation.json`, so a later rerun without a new evaluation does not change it. It is
   empty when a sample had no report (another filter route, or the reports were removed) and for evaluations published
   before the total was recorded; `samples.tsv` keeps the reads copied (`--reads`) or declared at preparation.

## Datasets

Declarations live in `datasets/`.

**Synthetic** (`kind: synthetic`): miRNAs with reads per sample (`mirnas`), miRNAs drawn from those annotated on the
reference (`random_mirnas`), isomiRs (`isomirs`: `3p_trim_1`, `3p_add_A`, `5p_trim_1`) and random reads; reads are
generated from a fixed seed with the kit's adapter (and `umi` random nucleotides at each end for kits that use them),
so the workflow's own filtering runs as usual. Not stored. This is the quick check on a reduced reference
(`synthetic_mmu_chr19`).

**Genome-wide synthetic** (`kind: synthetic`, `synthetic.mode: genome`, #412; `synthetic_hsa_genome`): reads of known
origin on the full human references, generated by `synthetic_genome.py` from the miRBase 23.1 mature, hairpin and GFF
files the declaration fills into `EXT_DIR/<data>/mirbase/23.1/` (`benchmark.py fill`, once).
- **miRNAs:** 500 annotated miRNAs. Accessions with the same mature sequence count as one, named by their first
  accession sorted, as the evaluation collapses them. Pairs one substitution apart are taken first (at 100:1 either
  way, or with the partner left out). Each is a family of its own: its other near-identical relatives are absent and
  related to no other chosen miRNA, so a read on a relative has one source; in large families such as miR-548 most
  pairs are therefore left out. Then matures at several loci, then the rest at random, both without near-identical
  relatives.
- **Abundances:** by rank, a Zipf curve (top 10 about 60% of the reads, as in plasma) and a low tail of 1–50 reads per
  million.
- **Mixed samples:** 6 + 6 in two groups, with about 1M reads each. log2 μ = log2 base + log2 fold change × group
  (sparse: 15% of the miRNAs, some near-identical and multi-locus) + loading × nuisance (a factor correlated 0.3 with the
  group, shared by miRNAs of one miRBase cluster). Counts are negative binomial with dispersion 0.05 + 2 / expected count.
- **Reads of a miRNA:** canonical (about 80%), 3' trimmed by 1–2 nt, 3' extended by 1–2 hairpin bases, a non-templated
  A or U, or the 5' end moved by one. An isomiR that would equal another mature, or a read of any of these classes that
  another miRNA's loci could give, is not generated; nor is an added A or U that any of the miRNA's hairpins has next.
- **Sequencing errors:** substitutions at 0.1% per base and 2% at the first base, with binned NextSeq qualities, low on
  the substituted bases; a few random reads.
- **Probe samples:** one canonical carrier set of miRNAs outside the truth, plus one read class of truth miRNAs (each
  isomiR class, or canonical reads with at least one error), two replicates each; and a carrier-only control.

`prepare` writes the reads, `truth.tsv` (per sample, class and read sequence, with its number of substitutions),
`samples.tsv` (each sample's kind, group, probe class and true factor scores) and `mirnas.tsv` (each miRNA's
abundance, fold change, nuisance loading, dispersion, cluster, loci, near-identical partners and probe classes). The
project uses the standard references (no `db.mirflowz.custom`). `--reads N` sets the mean reads of a mixed sample.

**Synthetic RNA pools** (`kind: pool`), from Giraldez et al. 2018 (*Nat Biotechnol* 36:746, doi:10.1038/nbt.4183),
one declaration per pool and kit, as a project has one kit:

| Declaration | Pool | Runs |
|---|---|---|
| `giraldez2018_equimolar_truseq` | 1,152 RNAs at equal amounts (GEO GSE94584) | 3 TruSeq runs, three labs |
| `giraldez2018_equimolar_nextflex` | same | 4 NEXTflex (4N) runs, one lab |
| `giraldez2018_ratiometric_truseq` | pools A and B, 334 RNAs at designed ratios from 1:10 to 10:1 (GEO GSE94585) | 2 A + 2 B per lab, two labs |
| `giraldez2018_ratiometric_nextflex` | same | 2 A + 2 B, one lab |

Each names the runs (ENA FASTQ with MD5), the composition table from the study's code
([rspengle/CrossU01_exRNA_Manuscript2017](https://github.com/rspengle/CrossU01_exRNA_Manuscript2017), MIT, at a fixed
commit) and the miRBase mature sequences, with SHA-256. `benchmark.py fill` downloads them once into
`EXT_DIR/<data>/` (each download to its own temporary file, renamed after the checksum matches, so fills of
declarations sharing a folder can run at the same time); it is the only step that writes to `EXT_DIR`:

```bash
python smallrna/benchmark/benchmark.py fill --dataset smallrna/benchmark/datasets/giraldez2018_equimolar_truseq.yaml \
    --ext-dir "$EXT_DIR"
```

`prepare` links the runs into the project, or copies the first `--reads N` reads of each for a quick run; `--methods`
picks a subset of the declared methods (mirflowz and TEsmall on human need the full human reference). Run `quant_all` with
`--notemp`: the oracle count reads the workflow's filtered reads, those of the project's filter route (`filter/cutadapt/`
with the kit default; the spike-in prefilter's when the project's `filter.spikein.ref` lists sets removed by Bowtie,
although the trimmed reads it starts from remain; the trimmed reads with native spike-ins, which UNITAS and miRge3
read), and the filter-route count reads those after each of the route's passes. Reads left by two trimming routes stop
`evaluate`.

**Another miRBase release** (#431). A declaration lists other releases under `mirbase_releases`, each with `url` and
`sha256` (and `hairpin` and `gff` for the genome-wide synthetic set); the human declarations list miRBase 22.
- `fill --mirbase-release 22` fetches that release's files into `EXT_DIR/<data>/mirbase/22/`.
- `prepare --mirbase-release 22` makes the project use it. The truth, the names and the release check come from it, and
  the project's `db.mirbase.release` is set to it for every method; `prepare` always sets the declaration's release.
- Every method's reference must then be built from that release:
  - UNITAS and miRge3 need a dump and a library built from it (`db.unitas.release`, `db.mirge3.release`). UNITAS stops at
    its reference check before counting; miRge3 records its library's release, and `evaluate` stops on a mismatch after
    quantification;
  - nf-core/smrnaseq and biocore_miraligner build theirs from gcfdb's miRBase files;
  - the genome routes need a genome reference for that release's assembly (#431).
- Evaluate each release in its own project, and compare across releases only where the miRNA set does not matter.

**Biological samples** (`kind: biological`): runs and metadata, no known composition, so the methods are compared with
each other, and agreement is not accuracy. A declaration lists `samples`, each with its `group` and `runs` (with `lane`
when known; `prepare` checks it against the read names), and the miRBase release for names and accessions. A sample's
runs are merged by the workflow (`R1` lists them comma-separated); `prepare --lanes` makes each run its own sample
(`<sample>_L<lane>`, so one run per lane and sample) for the lane-replicate metric. Run IDs must be unique and lanes
whole numbers; a lane is checked against the first 1,000 read names when they carry it.

| Declaration | Samples | Runs |
|---|---|---|
| `gse317995_nextflex_v4` | human plasma, NEXTflex v4 (GEO GSE317995): the two lowest-numbered samples of each of three groups | 4 per sample: the 4 lanes of one NextSeq 500 run |

GSE317995 is the only public NEXTflex v4 dataset found. Its lanes are one library each, so they replicate sequencing
only; library-preparation reproducibility needs separate libraries of one sample (e.g. the Giraldez 2018 plasma pool
made in several labs). `blockers: [blood miRNA]` is inferred, not recorded: miR-486-5p and miR-451a are nearly absent.

## Methods

`adapters.py` turns each method's native output into the common table: method, sample, the method's miRNA name, the
miRBase accession and the count (canonical reads and isomiRs together). Matching is by accession, because the methods
name miRNAs differently: UNITAS drops the species prefix and names reads it placed on a precursor rather than on a
mature miRNA (reads with a modification, e.g. a substituted first base or a shorter 3' end) after the hairpin and arm
(`miR-16-1-5p`), miRge3 (the workflow's own `miR.Counts.csv`, built from
miRge3's `mapped.csv` with the project-level exact-read rule, #349) reports merge groups (`name1/name2`, counted on the
first member), mirflowz names miRNAs after their locus (`mmu-miR-194-2-5p`) and reports isomiRs as separate rows,
TEsmall reports mature miRNAs by miRBase name (a name at several loci counted once), and tiny-count reports one row per
locus by miRBase feature ID (`MIMAT0000416_1`), summed per accession. Only TEsmall's `miRNA` rows are
scored: reads it assigns to the precursor (`hairpin`: reads lying less than 90% within the mature miRNA, e.g. a moved
5' end) are not mature counts and are left out. nf-core/smrnaseq's `mirna.tsv` names miRNAs by miRBase name, isomiRs
summed. biocore_htseq's `htseq_counts.tsv` (the genome counter of an earlier in-house pipeline, #425) names miRNAs by miRBase
name, isomiRs summed; a read counts in full for every miRNA name it overlaps, so miRNAs with identical sequences, scored
as one group, sum to a multiple of their reads. biocore_miraligner's `mirna.tsv` (#426) names miRNAs by miRBase name, its
substitution-free isomiRs summed, a sequence annotated to several miRNAs divided between them. A new method needs one
function in `ADAPTERS`, next to its smallrna
route.

A UNITAS name is placed by name when it is a mature miRNA of the organism, else by hairpin and arm: the hairpin of that
name in the UNITAS dump (`miRBase_hairpin_<Species>.fas`), and the organism's mature miRNA on that arm, the hairpin
family's (`miR-16-1-5p`: hsa-miR-16-5p) or the only one. Other names stay unplaced: no mature miRNA on that arm
(human miRBase has only hsa-miR-325-5p, so `miR-325-3p`), several of other families, or no such hairpin. Stripping
the locus number is not enough: hsa-miR-181a-2-3p is its own mature miRNA. On the Giraldez equimolar pool (TruSeq) this
places 4% of UNITAS's miRNA reads that had no accession and leaves 0.2% unplaced; the reads placed are those the other
methods count on the same miRNAs (often a pool sequence one or two bases shorter).

**isomiR table** (`isomirs.tsv`, #454). For the methods that report read sequences, one row per method, sample, reported
miRNA and read sequence, with the method's own label (canonical or isomiR, none for UNITAS). Each reader counts as the
method's own route does, so a miRNA's rows add up to its count in the common table; the evaluation stops when they do
not.
- **UNITAS:** `isomir_counts.tsv`, the name without the modification UNITAS appends to a read with a substitution, as
  `mir_counts.tsv` counts it.
- **miRge3:** `mapped.csv`, names folded by the project library's merge file, and only miRNAs with at least 2 exact
  reads in the project (as `miR.Counts.csv`); a read under `exact miRNA` is canonical.
- **mirflowz:** its isomiR rows. A plain row counts the reads matching mirflowz's annotated mature exactly; the table
  does not give their sequence, so they cannot be matched by sequence.
- **nf-core/smrnaseq:** `joined_samples_mirtop.tsv`, the first row of each UID; `Variant` NA is canonical.
- **biocore_miraligner:** each sample's `isomir.tsv`, the rows without a substitution, as its mature table.

TEsmall, tiny-count and biocore_htseq report features only and have no isomiR metrics. `isomirs.py` describes a
sequence relative to its mature miRNA (`variant`: 5' and 3' offsets up to 3 nt, up to 3 non-templated 3' bases, one
substitution; a templated extension is told from an added base only on the hairpin, so pools and biological samples
report `3p_longer`).

## Metrics

**Genome-wide synthetic data** (`genome_metrics.py`), per method, with `truth` (the truth counts analysed as a method,
the best the design allows) as a column:
- **Mixed samples:** `recovered`, `precision` (the smaller of truth and observed over the observed, so overcounting
  and reads on absent miRNAs lower it), `exact`, `misassigned` and `unplaced` (shares of the observed reads),
  `near_spill` (reads above the truth on near-identical partners, relative to the high member's truth) and
  `deviation_correlation` (each miRNA's log counts around its mean, observed against truth).
- **Differential expression:** median-of-ratios size factors, log2(normalised + 1), a linear model per miRNA with the
  group (and, `_adjusted`, the true nuisance scores), a moderated t (the residual variance shrunk towards the median,
  4 prior degrees of freedom), Benjamini-Hochberg at 5%. Each method tests the miRNAs with at least 10 reads on
  average in its own counts, as an analysis of them would (`de_tested`), so reads spilled onto a relative of a DE miRNA
  can make a false call. Reported as `de_sensitivity`, `de_fdr` (calls on any miRNA not designed DE, absent ones
  included; the report lists them) and `de_fc_error` (median |observed - true log2 fold change| of the true DE
  miRNAs). The ratiometric pools test edgeR and limma-voom (#414, Pools below).
- **Structure:** `structure_group_r2` and `structure_nuisance_r2` (each true factor regressed on the top two principal
  components) and `coexpression` (agreement of the miRNA-miRNA correlations).
- **Probe samples:** `probe_<class>`, the recovery of each read class, control subtracted. A generated read can still be
  a variant of another mature beyond these classes (a 2-nt 3′ trim of miR-320d-3p is the 5′ 17 nt of the miR-320
  family). The truth keeps its origin, which no counting method can know, so a method may count such reads on a
  relative.
- **isomiRs** (#454, mixed samples, no intervals yet): a true read counts as reported when the method reports its own
  sequence on its miRNA, up to the true count. `isomir_recovered_<class>` per read class (`canonical`, the isomiR
  classes, and `substituted` for reads with a substitution), and `isomir_recovered` over the isomiR classes. Over a
  method's reads: `isomir_misassigned` (on a sequence the sample has under another miRNA or as a random read) and
  `isomir_not_in_sample`. A method's own labels against the truth: `isomir_as_canonical` and `canonical_as_isomir`.

The report also gives recovery by stratum (unique, near-identical, multi-locus, identical sequences, low abundance),
each probe class's misplaced reads, and the miRNAs absent from the truth that receive reads.

**Synthetic data.** For truth as reads per miRNA, per method: `recovered` (share of the true reads counted on the right
miRNA), `exact` (share of miRNAs counted exactly, to within 1e-6 for the fractional counts of divided multimappers),
`misassigned` (reads on miRNAs not in the truth) and `unplaced` (reads on names without a miRBase accession).

**Pools** (`pool.py`):

- **Truth groups.** As in the study, sequences of 16-25 nt with only a 5' phosphate and unique in the pool are
  analysed, grouped against the organism's miRBase mature sequences: `exact` (the sequence of one or more mature
  miRNAs), `shifted` (a mature miRNA with up to 3 nt more or less at either end: expected as an isomiR), `foreign`
  (other species, viral, revised or non-miRNA sequences) and `excluded`. Against human miRBase 23.1 the equimolar pool
  has 407 `exact` sequences, the ratiometric pools 211.
- **Oracle.** The workflow's filtered reads counted by exact match to each pool sequence: what a sample contains after
  the kit's bias. Methods are scored against it on the `exact` group, so a kit's per-sequence bias (up to 100-fold) does
  not hide the differences between methods, and on the `shifted` group: isomiRs of known amount, which a method should
  count on their miRNA. An `exact` or `shifted` sequence whose miRNA is also that of another one is left out (in the
  equimolar pool 7 `exact` and 11 `shifted`, none in the ratiometric pools): a method's count of the miRNA holds the
  reads of both sequences, the oracle those of one. In `per_sequence.tsv` (oracle and method counts per `exact` and
  `shifted` sequence) such entries keep their oracle count, their method counts are `NA`, and `shared_with` names the
  others; `shift` gives the kind of shift.
- **Kinds of shift** against the miRNA's mature sequence: `3p_shorter` and `3p_longer` (the 5' end and so the seed in
  place), `5p_moved` (another seed, whatever the 3' end does) and `mixed` (a sequence shifted from several miRNAs in
  different ways). The report gives the `shifted` group's metrics per kind.
- **Equimolar pool**, per method and for the oracle: `detected`, `cv` and `qcd` (quartile coefficient of dispersion)
  across the group; per method against the oracle: `median_log2_vs_oracle` (of (method + 1) / (oracle + 1)),
  `within_2fold`, `beyond_10fold`, `spearman_vs_oracle`, and of all the reads the method reported as miRNAs, `outside`
  (the share on miRNAs of no `exact` or `shifted` sequence, i.e. assigned from foreign sequences) and `unplaced` (the
  share on names without a miRBase accession). Over the `shifted` group, with a `shifted_` prefix: `detected` and the
  metrics against the oracle.
- **Ratiometric pools**: counts per million of the group, averaged per pool, give an observed log2(A/B) per sequence,
  centred on the sequences designed 1:1. Against the design: `median_abs_error`, `spearman_vs_design`, `called`
  (designed at least 2-fold apart and observed so, in the right direction) and `false_calls` (designed less than 2-fold
  apart, observed at least 2-fold apart). A fold change cancels the per-sequence kit bias, so what differs between
  methods is quantification. The report also lists the median observed log2(A/B) per designed ratio. Over the
  `shifted` group, with a `shifted_` prefix, the same four, its sequences scaled by the `exact` group's counts and
  centred on its 1:1 sequences, so the `exact` results do not depend on them.
- **Differential expression on the ratiometric pools** (#414; `de:` and `evaluate --de`). Each method's counts of the
  scored `exact` sequences, rounded to whole numbers, are tested pool A against pool B with edgeR's exact test, edgeR's
  GLM likelihood-ratio test and limma-voom. Library sizes are all the miRNA reads the method reported in a sample, the
  workflow's counts per million (the oracle: its pool reads). There is no TMM and, unlike the metrics above, no
  centring, so an offset between the pools stays in the fold changes. The exact test runs as the best configuration in
  the McCall lab's benchmark of miRNA-seq normalisation and differential expression (Aparicio-Puerta et al., bioRxiv
  2026, doi:10.64898/2026.05.08.723923): a common dispersion and the exact test, so two groups without covariates. The
  GLM test and voom take the declaration's covariates (`de: {covariates: [lab]}`: `~ lab + group`). The GLM test keeps
  the exact test's single common dispersion but estimates it under that design. Without covariates it calls nearly what
  the exact test calls, so it differs from the exact test mainly by the covariates, and from voom by the method. DESeq2
  is left out: it did worst there. A call is a BH-adjusted p of at most 0.05. Metrics, prefixed `edger_`, `edger_glm_`
  and `voom_`: `called` (sequences designed at least 2-fold apart called in the designed direction), `false_calls` (1:1
  sequences called; a call on a sequence designed 1.5-fold apart is a true one), `fdr` (the share of calls on 1:1
  sequences or in the wrong direction), `fc_bias` (median of sign(d) × (log2 fold change − d), negative when fold
  changes are compressed) and `offset` (median log2 fold change of the 1:1 sequences). Against the oracle's fits on
  the same runs, for each method (#444): `own_false_calls` (1:1 sequences the method calls and the oracle does not),
  `discordant_vs_oracle` (sequences whose call, up, down or none, differs from the oracle's), and over the sequences
  both test `fc_shift_vs_oracle` (the median difference of the log2 fold changes, method − oracle) and
  `fc_error_vs_oracle` (their median absolute deviation from that shift). The shift is common to all sequences and
  comes from the library sizes, all of a method's miRNA reads against the oracle's pool reads; it also moves calls near
  the cutoff. The error is the disagreement of the counts themselves. The pools' 1:1 sequences are not all 1:1 in the
  libraries: the oracle calls some of them, and different ones in each lab of the TruSeq pools. So `false_calls`
  measures the libraries as well as the method, while the metrics against the oracle compare the methods on the same
  libraries.
  Their intervals resample the sequences on the fixed fits, on every draw: a draw without a 1:1 sequence has no
  `offset` or `false_calls` but still scores the rest, while the metrics above skip it. The declaration must leave the
  design residual degrees of freedom. The report gives the share called per designed fold change and the tools'
  versions and settings; `de/de.log` has each row's two common dispersions.

  With `de: {within: lab}` (#444) each lab's runs are also tested on their own with `~ group`, as replicate benchmarks
  of a single-lab experiment; each lab needs at least two runs per pool. Their metrics carry the lab, e.g.
  `edger_fc_error_vs_oracle.Lab1`. They are in `metrics.tsv`, `metrics_ci.tsv` and `differences.tsv`, whose paired
  differences show whether a difference between methods holds in each lab. The report shows their values in a section
  of their own, outside the main table and the uncertainty section.
- **isomiRs** (#454, no intervals yet). Against the oracle, per run: `canonical_isomir_recovered` and
  `shifted_isomir_recovered`, the share of the oracle's reads of the scored `exact` and `shifted` sequences a method
  reports as that sequence on the entry's miRNA, and `shifted_isomir_within_2fold`. A synthetic pool holds only its
  designed sequences, so a method's reads on other sequences of its miRNAs are library or sequencing artefacts
  (`artefact_share`, broken down by kind in the report); all methods read the same filtered reads, so the differences
  are the methods'.
- **Filter route.** The same exact-match count of the pool sequences after each pass of the project's filter route
  (`filter_stages.tsv`): for cutadapt after the 3' adapter, the filter (random ends, 5' adapter, length, N, expected
  errors) and the artefact pass, for fastp after the adapter and the filter pass, then the spike-in prefilter with
  Bowtie spike-ins. The reads after the adapter pass still carry the kit's random ends; they are cut before matching, by
  as many bases as the next pass cuts (cutadapt `-u`, fastp `-f`/`-t`, read from its JSON report). The report gives the
  share of the pool reads each later pass removes, per truth group, and the sequences with at least 20 reads after the
  adapter pass that lose more than 10% of them, with their length and GC content. A read that is not an exact pool
  sequence after the adapter pass (a sequencing error, an untrimmed adapter) is counted in no pass; the reads before the
  adapter pass are not compared.

**Biological samples** (`biological.py`; detected = at least 10 reads in a sample; CPM relative to each method's
miRNA reads):

- **Accessions.** miRNAs with the same mature sequence count as one (the lowest accession), since the methods credit
  their reads to different members (miRge3 to its merge group's first name; UNITAS and mirflowz to either, e.g.
  hsa-miR-199a-3p and hsa-miR-199b-3p). Reads on names without a miRBase accession are left out of every metric and
  reported as `unplaced`, with the names that have most.
- **Per method:** miRNA reads, `unplaced`, miRNAs detected per sample and in every sample, and the share of the 10 most abundant.
- **isomiRs** (#454), a description: `isomir_share` (the share of a method's miRNA reads that are not the mature
  sequence) and `isomirs_detected` (sequences with at least 10 reads over all samples); the report adds each method's
  reads by kind and the overlap of detected isomiRs between methods (Jaccard).
- **Between methods**, per pair (medians over samples): Spearman of log2(CPM + 1) over the miRNAs either detects
  and over those both detect, median |log2 CPM ratio| over those both detect, overlap of the 100 most abundant, and miRNAs only one detects
  (`only_a` and `only_b` for the first and second method of the pair).
- **Blocked miRNAs:** CPM of miR-451a, miR-486-5p and miR-92a-3p, the targets of the NEXTflex blood miRNA blockers.
- **Filtering:** reads through the cutadapt route's three passes (3' adapter found, 5' adapter discarded, too short,
  N, expected errors, 3' artefacts and reagent oligo reads discarded, kept), from their JSON reports, which must chain (each pass read what the
  one before wrote); only for the cutadapt route.
- **First bases:** in the 20 most abundant miRNAs, the share of reads that match from the second base on but have
  another first base, with the median quality of canonical and substituted first bases (reads of miRNAs that differ
  only at the first base, such as hsa-miR-17-5p and hsa-miR-106a-5p, keep their own miRNA; a third first base is not
  placed). On NEXTflex v4 (no random 5' bases, so the insert starts at cycle 1 of a two-colour NextSeq) about 2-3% have a
  substituted first base of low quality: a sequencing error, in every lane but at a rate that differs by lane.
- **Lane replicates** (`--lanes`): per method, the dispersion index of lane counts against sampling alone
  (multinomial within the method's miRNA reads; miRNAs with at least 20 reads per lane), with the values sampling
  alone gives (mean 1; with 4 lanes a median of 0.79 and about 3% above 3), how many miRNA–sample pairs exceed 3
  (`mirna_sample_pairs` tested), that share by mean reads per lane, and the bin of mean reads per lane from which the
  lanes' CPM agree within 10% (median sample CV). All methods count the same reads, so a lane effect raises every
  method's values, and most for abundant miRNAs, whose sampling noise is small: a share above 3 that grows with
  abundance is a multiplicative lane effect (e.g. an error rate that differs by lane), not a method's. Compare the
  methods with each other rather than with the sampling values.

On the mouse chr19 dataset, UNITAS and miRge3 count every miRNA exactly. mirflowz counts a third of
`mmu-miR-146b-5p`, whose sequence occurs three times on chr19: it divides a read that maps equally well to several loci
between them, while UNITAS and miRge3 match by sequence.
