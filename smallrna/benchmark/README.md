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
   With `--mirgenedb GFF3` (#493; MirGeneDB's annotation of the organism, e.g. miRge3's library file
   `annotation.Libs/<species>_MirGeneDB.gff3`), every miRBase accession gets a stratum: `mirgenedb` (listed), `star`
   (listed as a star arm) or `mirbase_only`. Identical sequences count as one, in the first stratum any of them has.
   - **Biological datasets:** `detected.<stratum>`, `detected_alone.<stratum>` (detections no other route makes in
     the sample) and `reads.<stratum>` are added to the metrics, and the report lists the miRNAs the routes dispute
     most (totals at least 4-fold apart), with their stratum.
   - **Checks:** the file's MirGeneDB IDs must be the dataset organism's (`Hsa-`, `Mmu-`, `Rno-`), and the report names
     the release from its header.
   - **Null comparisons (#449):** `edger_null_called.<stratum>` and `voom_null_called.<stratum>`, the null splits'
     false calls among each stratum's tested miRNAs.
   - **The mixture series:** `linearity_error.<stratum>` and `within_1_5fold.<stratum>`.
   - **The genome-wide synthetic set:** the strata become rows of its recovery-by-stratum table.
   Before scoring, it sets a method aside as `stopped` when its count table is missing or does not have exactly the
   project's samples (a missing sample would count as zero), or, for pools and biological datasets, when the method's
   `quant/<method>/reference.tsv` is missing or names another miRNA database or release than the dataset's
   `mirbase.release` (#407); the other methods are scored without it (#448). It stops when the project's configured
   samples differ from the benchmark's, or when no method can be scored. UNITAS and miRge3 record their
   reference's release; `quant_reference` records it for mirflowz, TEsmall, tiny-count and biocore_htseq (`db.mirbase.release`,
   or `NA` for a `db.mirflowz.custom` reference) and nf-core/smrnaseq and biocore_miraligner (`db.mirbase.release`), with the
   SHA-256 of the count
   table it was written after, which must match the table.
   **RNA classes and read lengths** (#504, classes.py), for every dataset kind: each evaluated route's own class
   assignment, read from the files it writes per read sequence, mapped to common classes (miRNA, miRNA of other
   species, hairpin, tRNA, rRNA, snoRNA, snRNA, misc ncRNA, piRNA, lncRNA, mRNA, intron, pseudogene, TE, low
   complexity, other, `unmapped` for a native class the mapping does not know, `unassigned`). Nothing is quantified
   again: UNITAS's `unitas.full_annotation_matrix.txt` (a sequence divided evenly between its annotations, which
   reproduces its `annotation_summary.txt`), miRge3's `mapped.csv`, TEsmall's `.anno.rlen.info`, `.rRNA.bam` (the reads
   its rRNA filter removed) and `count_summary.txt` (structural RNA by type, totals only), tiny-count's
   `assigned_nt_len_dist.csv`, smrnaseq's and miraligner's isomiR tables, and biocore_htseq's `isomir.tsv` and, for
   human, `ncrna.tsv` (#503). The baseline is each sample's filtered reads by length, from miRge3's `mapped.csv` and
   `unmapped.csv`, which hold every filtered read, or else the filtered FASTQ; `unassigned` is the baseline minus the
   reads a route assigned (biocore_htseq: the sequences in neither of its tables). Per route and class
   `class_states.tsv` gives `reported`, `lumped` (miRge3's `ncrna others`, TEsmall's structural RNA and gene exons),
   `totals only`, `not reported` or `not available` (a file the route writes is missing, e.g. TEsmall's rRNA BAM:
   its rRNA and unassigned reads are then not given); the report shows the median share of each sample's filtered reads per class and
   route, every row by sample and length is in `classes.tsv`. Where matplotlib is installed, two plots follow MultiQC's UNITAS section, pooled over samples as % of
   the filtered reads, with the classes grouped (miRNA, tRNA, rRNA, snoRNA, misc ncRNA with snRNA, lncRNA, mRNA, other,
   lumped, unassigned): `classes.png`, one stacked bar per route, and `class_lengths.png`, reads by length per class,
   one line per route that reports the class. Classes come from each route's own references and precedence, so their shares
   partly reflect those. A route whose files cannot be read is left out with a note; the evaluation goes on.
   On a synthetic set with decoy probe samples (#412), `class_truth.tsv` and a report table score the class
   assignment against the truth: each `P_decoy_<class>-<source>` sample adds the reads of one class to `P_control`'s
   carriers, and the share of those reads (that passed the filter) a route puts in the right class, in a coarser
   class holding it, in miRNA, elsewhere or nowhere. Routes with outputs per sequence are read again without the
   truth's carrier sequences (a route can place a carrier differently from sample to sample); TEsmall and tiny-count
   by their classes minus `P_control`'s. Right classes: tRNA for GtRNAdb and mitochondrial tRNA, rRNA, misc ncRNA for
   Y RNA, snoRNA, snRNA or misc ncRNA for `other_ncRNA` (snRNA, 7SK, 7SL, vault RNA), mRNA, piRNA, and hairpin for
   miRNA precursors. The rRNA decoys are also scored by molecule from the truth's origin, for the routes with outputs per
   sequence: RefSeq's 18S, 5.8S, 28S, 5S and the 45S pre-rRNA's spacers (ETS, ITS; rare in small RNA libraries),
   Ensembl's 5S, 5.8S and mitochondrial rRNA. The decoys end at the annotated boundaries, so fragments such as the 5.8S
   long form's 5' end are not in the truth.
   Every metric gets a 95% bootstrap interval (#413), and so does every difference between two methods, paired
   because both are computed on the same resample. A resample draws the scoring units (truth miRNAs; the scored exact
   and shifted pool sequences) and, where there are several, the runs or samples, the ratiometric pools' runs within
   pool A and pool B, with replacement. `metrics_ci.tsv` (metric, row, value, low, high) and `differences.tsv`
   (metric, row_a, row_b, difference, low, high) are published with the metrics and listed with their SHA-256 in
   `evaluation.json`, with the resamples and seed (`--bootstrap N`, default 1000, `0` for none; `--seed`, default 1).
   `report.md` has a section with the intervals and the differences whose interval excludes 0. Some values have no
   interval: the synthetic totals `misassigned` and `unplaced`, and biological `detected_all_samples`, which a
   resample repeating a sample biases upwards. With few samples (e.g. the 6 of the GSE317995 quick run) the sample part of the interval
   is approximate.
   **Route states** (#448): every declared route is `evaluated`, `not run` (not selected with `prepare --methods`,
   recorded as `not_selected` in the project's `benchmark/dataset.yaml`; or a run not done, e.g. for resources,
   recorded by hand in `benchmark/not_run.tsv` with the columns `route` and `reason`) or `stopped` (a check above
   failed, with the check's message). `routes.tsv` and the report's Routes section list them, and the summary shows
   the state in the route's cells.
   **Counting settings** (#417): each route records the settings that decide its counting with its count table
   (`quant/<route>/settings.tsv`, written by the workflow's `quant_settings` rule from the values its commands use, with
   the table's SHA-256): e.g. tiny-count's mismatches and hit normalisations, TEsmall's class order, miRge3's
   arguments and project rule, mirflowz's run configuration. The report's Routes section lists them per evaluated
   route (`fixed` ones set in the rule or script), and `settings.tsv` holds them all. A record written for another
   table stops the route; a project quantified before #417 has none ("not recorded").
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

   **Differential expression** (#414, #447). One step for every dataset that declares `de:`: a ratiometric pool, the
   genome-wide synthetic dataset or a biological dataset. `evaluate` also writes the inputs to `benchmark/de/input/`:
   the counts of each row (a reference row, the pools' oracle or the synthetic truth, and each method), the library
   sizes, the design and the model (`model.tsv`: the contrast, the covariates that are numbers and the least mean
   count a row tests). A dataset can have more DE sets, each in a folder of its own, listed in
   `benchmark/de/levels.tsv`:
   - a pool or biological dataset with `de: {within: lab}` (#444): each lab's runs, in `benchmark/de/lab/<lab>/`;
   - the synthetic dataset: the same samples with the true nuisance scores as a covariate, in
     `benchmark/de/adjusted/nuisance/`.

   `de.smk` tests every set with edgeR's exact and GLM tests and limma-voom (`de.R`, in the `sartools` image), and
   `evaluate --de` adds the scores of their calls (Metrics, below):

   ```bash
   snakemake -s smallrna/benchmark/de.smk -d <project> --use-apptainer --cores 1
   python smallrna/benchmark/benchmark.py evaluate --project <project> --de
   ```

   Run `de.smk` with the same container settings as the project's `quant_all`: the repository holding `de.R` must be
   visible in the container. `de.R` fits copies of the inputs and records their MD5s, its own and that of the results
   table it writes. `evaluate --de` stops before writing anything unless that record matches the current inputs, `de.R`
   and the results table; after new counts, run `evaluate` without `--de` (it writes the current inputs), then
   `de.smk`, which also reruns after a change to `de.R`. A project prepared before its dataset declared `de:` enables
   DE by adding the dataset's `de:` to `benchmark/dataset.yaml`. The inputs changed with #447 (`model.tsv`, and the
   column `feature` in the counts and results), so a project evaluated earlier needs `evaluate` once more before
   `de.smk`. The DE files are not listed in `evaluation.json`; the report gives each set's results SHA-256.

   What a declaration's `de:` takes:

   | Dataset | `de:` | Contrast | Rows test |
   |---|---|---|---|
   | ratiometric pool | `covariates`, `within` | pool A against pool B | every scored sequence with reads |
   | genome-wide synthetic | nothing (`de: {}`) | the group the designed fold changes raise against the other | miRNAs with at least 10 reads on average |
   | biological | `contrast: [numerator, denominator]`, `covariates`, `within`, `null` (#449) | two of the samples' groups; samples of other groups are left out | miRNAs with at least 10 reads on average |
   | mixture | nothing (`de: {}`) | the pure samples against each other, and each pair of the mixtures' fractions, the higher against the lower; each paired by subject | miRNAs with at least 10 reads on average |

   A covariate is any field of the declaration's runs (pools) or samples (biological), e.g. `lab`, `patient` or
   `mouse`; each needs a value for every run or sample of the contrast, and the design must leave residual degrees of
   freedom. Group labels and covariate values are text, whatever YAML makes of them (a dose `100` against `0`). A row
   with a single feature to test is not tested: one feature gives no common dispersion. A biological dataset prepared
   with `--lanes` has no DE: the lanes of one library are no replicates. Its calls are reported, not scored, since no
   change is known; with `within`, the report has a table for the whole contrast and one for each level.

   **Null comparisons** (#449, `de: {'null': {splits: 50, seed: 1}}`; the key quoted, since a bare `null` is YAML's null): the samples of one group split in two at random
   and tested against each other, so every call is false: the false-call rate under real variance. Each group with at
   least 4 samples gets up to `splits` distinct splits drawn from the seed (a split and its mirror are one), each a DE
   set (`benchmark/de/null/<group>-<k>/`) with pseudo-groups `x` and `y`. A split halves each stratum of the declared
   covariates and `within`, so the halves are balanced on the design; a covariate the group's samples share is left out
   of its splits. Per method and test, over all splits: `<test>_null_called` (the median share of tested miRNAs
   called, the report adds its range), `<test>_null_any` (the share of splits with a call) and `<test>_null_p05` (the
   share of tested p-values of at most 0.05, 5% when they are uniform); the report also lists the miRNAs called in the
   most splits. Their spread over the splits is their uncertainty, so the sample bootstrap leaves them out.
5. **Summary.** `benchmark.py summary` reads several evaluated projects and writes one overview, a Markdown file or a
   single self-contained HTML file (no external resources), by the output's extension:

   ```bash
   python smallrna/benchmark/benchmark.py summary --projects <project> [<project> ...] --output summary.html
   ```

   It has a section per dataset kind (synthetic, equimolar pool, ratiometric pool, biological) with a row per project:
   its dataset, kit, samples, reads processed, methods, a link to its
   `report.md`, and the kind's key metrics from its `metrics.tsv` per method (and the oracle, where it has values),
   each with its bootstrap interval in parentheses when the evaluation has one.
   A declared route that was not evaluated shows its state, `not run` or `stopped` (#448); a metric missing from a
   project's table, or a method not declared there, is an empty cell.
   Projects of one synthetic declaration that differ only by the seed of the draw (`prepare --synthetic-seed`, #455)
   also get a table over the seeds after their section: per key metric and route the mean and the range, and how many
   pairs of the routes with a value in every seed keep their order in every seed. Two projects with one seed stop it. The seed is read from `evaluation.json`
   (`synthetic_seed`), or from the project's `benchmark/dataset.yaml` for an earlier evaluation.
   A settings variant (`prepare --set`, #417) gets a table against the project of the same dataset, samples, draw and
   reads processed run with the defaults (metrics over the routes' shared features left out): per key metric, the
   varied route's value with the defaults and with the setting, and the routes whose order with it (higher, lower or
   equal) the setting changes. Two such default projects stop it. It writes nothing into the
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
6. **Depth.** `benchmark.py depth` (#409) gives the sequencing-depth (saturation) curves of an evaluated biological
   or mixture project, from its tables, without running the workflow again:

   ```bash
   python smallrna/benchmark/benchmark.py depth --project <project> [--output <dir>] [--depths 500k,1M,2M,4M,8M,16M,32M] [--seed 1] [--plot]
   ```

   Each method's common table (and isomiR output) is thinned binomially to each depth: every read is kept with
   probability depth / the sample's sequenced reads (the first cutadapt pass's input, from the reports kept with
   `--notemp`), the same draw for every method. For detection this equals subsampling the reads when a read counts
   for one miRNA; tiny-count's and mirflowz's fractional counts are rounded first, so it is approximate for them.
   Depths above the shallowest sample are left out, and each sample's full depth ends the series. Merged samples only.
   - **Per method and depth** (`depth.tsv`, `depth.md`): median miRNA reads per sample; miRNAs with at least 1, 5 and
     10 reads (median per sample) and with at least 10 in every sample; the replicate CV, the median over groups of
     the median CPM coefficient of variation of the miRNAs with at least 10 reads in every sample of the group; and
     the distinct read sequences of the isomiR output (median per sample).
   - **Per pair of methods** (`depth_agreement.tsv`): the biological report's agreement (Spearman over miRNAs either
     detects, median |log2 ratio|, top-100 overlap) at each depth.
   - **Saturation:** the depth at which detection (at least 10 reads) reaches 90% and 95% of its full-depth value,
     interpolated on log depth.
   - **Plot** (`--plot`, needs matplotlib): `depth.png`, detection and replicate CV against depth.
   - **Differential expression by depth** (#414), for a dataset that declares `de:`. The evaluation's DE sets (the
     contrast, each level of `within`, the null comparisons of #449) are tested on each depth's thinned tables (the
     draw of the detection metrics) and on the full depth. The tests are de.R's, over the full depth's miRNAs, with
     every miRNA read a route reported at that depth as its library size. `depth` writes the inputs to `de/<depth>/<set>/`
     (`all` for the contrast, `lab.L1` for a level, `null.case-01` for a split), with the full depth's equal to the
     evaluation's. de.smk fits them, and `depth --de`, with the same depths and seed, checks and scores the results:

     ```bash
     snakemake -s <repo>/smallrna/benchmark/de.smk -d <depth output> --config de=de --use-apptainer --cores 8
     python smallrna/benchmark/benchmark.py depth --project <project> --output <depth output> --de
     ```

     Per route, test and depth (`depth_de.tsv`, `depth.md`), against the route's own full-depth calls (BH-adjusted p
     of at most 0.05, in the same direction): `calls`, `recall` (share of the full-depth calls made) and `precision`
     (share of the calls the full depth makes). For the null comparisons, `null_called`, `null_any` and `null_p05` per
     depth, as in the evaluation. Nothing is written when a result comes from other inputs (other depths or seed).
     These scores show how depth changes a route's calls, not whether the calls are right. They are optimistic: each
     depth is a subset of the full-depth reads, and the full depth is what the project holds (a read cap where one was
     set). de.smk fits each set at each depth, plus the full depth, one fit of a few seconds each. 50 null splits in
     each of three groups give 150 sets, so 600 fits over three depths and the full depth: run de.smk with `--cores`
     to match.
   - **A mixture** (#451, #488): its own DE sets and library sizes (reads on miRNAs with an accession). The pure
     samples are scored as above, against the full depth's calls. The contrasts between fractions are scored against
     each route's expected changes, which come from its unthinned pure samples and content ratios and are the same at
     every depth, with the evaluation's metrics: `called`, `called_shared`, `wrong_direction`, `fc_bias` and
     `unchanged_called`, set `between`. `scored` counts the miRNAs a route expects at least 2-fold apart and still
     tests, since `called` is a share of those.

   It writes to `--output`, or `benchmark/depth/` in the project; it changes nothing else there.
7. **Depth imbalance.** `benchmark.py imbalance` (#488) is a stress test of the normalisation arms for an evaluated
   biological or mixture project that declares `de:`, without running the workflow again. It thins each DE set's samples
   unevenly from the count tables, binomially as `depth` does, and tests the sets again with de.R:

   ```bash
   python smallrna/benchmark/benchmark.py imbalance --project <project> [--output <dir>] [--ratios 1.5,4] [--seed 1]
   snakemake -s <repo>/smallrna/benchmark/de.smk -d <output> --config de=de --use-apptainer --cores 16
   python smallrna/benchmark/benchmark.py imbalance --project <project> [--output <dir>] [--ratios 1.5,4] [--seed 1] --de
   ```

   - **Designs** per set (the contrast, each level of `within`, the #449 null splits), at each ratio r:
     - `independent-<r>`: half the samples of each group, a seeded draw, thinned to 1/r.
     - `tied-<r>`: the set's second group only (the contrast's denominator, a split's `y`) thinned to 1/r, so a
       contrast is tested with the imbalance in one direction.
     - `none`: unthinned.
   - **Arms:**
     - `workflow`: every miRNA read the route reported, the evaluation's library size.
     - `fixed`: the sample's filtered reads times its thinning fraction.
     - `tmm`: edgeR's TMM factors on the workflow's.

     Inputs go to `de/<design>_<arm>/<set>/`. `none` in the `workflow` arm equals the evaluation's inputs.
   - **Scores** (`imbalance.tsv`, `imbalance.md`), per design, arm, route and test:
     - the null splits' `null_called`, `null_any` and `null_p05`, and `offset`, the median over the splits of the
       median log2 fold change (every change is null, so any shift comes from the normalisation);
     - the contrast's sets against the same arm's unthinned fits: `calls`, `recall`, `precision` and `fc_shift`, the
       median change of the log2 fold changes;
     - a mixture's contrasts between fractions against the expected changes, as for `depth` (set `between`). There
       every miRNA has an expected direction and size, so this is the best real test of what an arm does. The expected
       changes are on the scale of a route's miRNA reads, the `workflow` arm's. `fixed` is scored against them shifted,
       per subject before the mean, by the log2 ratio of the route's miRNA share of filtered reads between the two
       fractions, from the unthinned data (the mean over the subjects is printed per contrast). `tmm` is scored on the
       workflow's scale, as TMM meeting a design where most miRNAs change. A mixture's `tied` thins each contrast's lower fraction, the samples with less of the second tissue.
   - **Realised depths** (`imbalance_depths.tsv`, a table in `imbalance.md`): per design and set, the median
     filtered reads times the fraction of each group, and their ratio. The ratios act on the project's own depths,
     which may differ between the groups already in `none` (a proportional read cap keeps them).
   - **Caveat:** thinning keeps each sample's composition and every arm's library sizes follow it. The designs
     therefore show the cost of the imbalance (power, false calls from the shallower samples' variance), not a bias
     from composition.
   - **Fits:** sets × (1 + 2 × ratios) × 3 arms. 150 null splits at two ratios make 2,250 fits.

   Nothing is written when a result comes from other inputs (other ratios or seed), and new inputs remove an earlier
   `imbalance.tsv`.

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
- **Decoys** (#412 PR 2, `synthetic.decoys`, decoys.py): reads of other RNA classes, from a random stream of their own
  after the miRNA reads, so those stay the same with or without decoys. Each class and source is a unit:
  - **rRNA:** RefSeq 18S, 28S, 5.8S, 5S and the 45S spacers; Ensembl mitochondrial 12S/16S and 5S/5.8S genes.
  - **tRNA:** GtRNAdb mature tRNAs, and Ensembl mitochondrial tRNAs, both with CCA added. Fragments are 5′ and 3′
    halves, 5′- and 3′-tRFs and i-tRFs, some with an m1G9 or m1A58 mismatch.
  - **Ensembl:** Y RNA, snoRNA and scaRNA, other ncRNA (snRNA, 7SK, 7SL, vault RNA), and mRNA.
  - **piRNA:** piRNAclusterDB clusters, 24–32 nt with a 5′ U bias.
  - **Precursor:** miRBase hairpin loops and flanks outside the matures.

  The sources are filled with checksums into `<data>/decoys/`.

**Scenarios** (#453): a declaration can extend another of the same folder (`extends: <file>`) and change some of its
values (mappings are merged key by key, anything else replaced; one level only, with a name of its own). A scenario
keeps only its changes, so a later edit of its base (a new route, another release) applies to it too; it may change
the design (`synthetic.structure`, `de`, `methods`), not the dataset or its references (`kind`, `organism`, `libprepkit`,
`data`, `mirbase`, `synthetic.mode`, `synthetic.decoys`). The generator draws the group design, the isomiR design and
read errors, and the probe samples from streams of their own, so a scenario keeps its base set's miRNAs, abundances,
nuisance loadings, isomiRs and probe samples at the same seed; only the mixed samples' counts follow its design. The
genome-wide synthetic set has two, with the same seed, data folder, decoys and probe samples:
- `synthetic_hsa_genome_null`: no miRNA designed to change (`structure.group.fraction: 0`), so every DE call is false:
  each route's, test's and normalisation arm's false-call rate under the generator's variance alone;
- `synthetic_hsa_genome_unbalanced`: every change goes up (1.5-, 2- or 4-fold), and the 5 most abundant miRNAs (about
  half the miRNA reads) always change (`structure.group.abundant: 5`, which makes the most abundant miRNAs change
  first; without it the draws are the base set's). The composition shifts, so library sizes made of all miRNA reads move
  the unchanged miRNAs the other way; `false_calls_against` counts the calls that makes, and the normalisation arms
  show how much `fixed` and `tmm` remove.
  - **Fragments:** each unit draws a fixed catalogue of fragments with lognormal weights.
  - **No look-alikes:** a fragment sharing 18 nt with a miRBase hairpin is drawn again, so these decoys measure reads
    unrelated to miRBase.
  - **Mixed samples:** `load` gives decoy reads per miRNA read, rising with the nuisance factor (`nuisance`).
  - **Probe samples:** each unit has its own (`P_decoy_<unit>_<n>`, kind `decoy`), with the carrier.
  - **miRBase matures inside other-class molecules** (e.g. tRNA- and rRNA-derived entries) are probed from both sides:
    their canonical reads (`P_embedded_<class>_<n>`), and fragments of the host molecules that cover at least half of
    such a mature without being a read any miRNA could give (`P_overlap_<class>_<n>`, kind `overlap`).

`prepare` writes:
- the reads;
- `truth.tsv`: per sample, class and read sequence, with its number of substitutions, its `source` and, for decoys,
  its `origin` (molecule and fragment type);
- `samples.tsv`: each sample's kind, group, probe class and true factor scores;
- `mirnas.tsv`: each miRNA's abundance, fold change, nuisance loading, dispersion, cluster, loci, near-identical
  partners, probe classes, and the other-class molecules it lies `inside`.

The project uses the standard references (no `db.mirflowz.custom`). `--reads N` sets the mean reads of a mixed sample.

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

`prepare` links the runs into the project, or copies the first `--reads N` reads of each for a quick run, or with
`--fraction F` the first F of each run's declared reads: the same total can be had either way, but `--reads` gives every
run the same depth, while `--fraction` keeps the depths the runs have relative to each other, which the normalisation
in differential expression is there to handle. `--methods`
picks a subset of the declared methods (mirflowz and TEsmall on human need the full human reference), and, for a
biological dataset, `--per-group N` the first N samples of each group in declared order; for a synthetic dataset,
`--synthetic-seed N` another draw of the same declaration (#455: the miRNAs chosen, their abundances and changes, the
counts and, for the genome-wide set, the decoys), each seed a project of its own, recorded in its
`benchmark/dataset.yaml` and in `evaluation.json`, and for the genome-wide set `--no-probes` the mixed samples alone,
without the probe, decoy, overlap and control samples (`no_probes: true` in the project's `benchmark/dataset.yaml`).
The mixed samples' reads and truth are the ones the same seed gives with them. Their metrics are the same for a route
that counts each sample on its own; a route with a rule over the whole project can differ on rare miRNAs (miRge3
keeps a miRNA with two exact reads anywhere in the project, and a probe sample can supply one of them). The project
has no `probe_*`, decoy or overlap metrics, its report says so, and the summary still groups it with the other seeds
of its declaration. It is
for further seeds and scenarios, where the probe samples (most of a project's samples, a few percent of its reads)
would test the same read handling again. `--set quant.ROUTE.SETTING=VALUE` (repeatable) a settings variant
(#417): a route setting of `smallrna.config` changed for the project, for a route it runs, not to its default. The
variant is recorded as `variant` in `benchmark/dataset.yaml` and named in the report; run it with only the varied route
(`--methods`) and compare it with the project run with the defaults in `summary`. #417's variants: tiny-count
`mismatch` 1 and `normalize_by_genomic_hits` false, TEsmall's own class `order`
(`structural_RNA,miRNA,hairpin,exon,TE,intron,piRNA_cluster`). miRge3's `--crThreshold` is no variant: the scored
table is rebuilt from miRge3's alignments (exact and isomiR reads of miRNAs with 2 exact reads in the project), which
the threshold does not reach; mirflowz's multimapper handling is not configurable. `evaluate` scores a variant only when
its route's settings record (`quant/<route>/settings.tsv`) has the setting at the variant's value. Every declaration lists the
eight routes (#448); the ones left out at `prepare` are `not run` in the evaluation. Run `quant_all` with
`--notemp`: the oracle count reads the workflow's filtered reads, those of the project's filter route (`filter/cutadapt/`
with the kit default; the spike-in prefilter's when the project's `filter.spikein.ref` lists sets removed by Bowtie,
although the trimmed reads it starts from remain; the trimmed reads with native spike-ins, which UNITAS and miRge3
read), and the filter-route count reads those after each of the route's passes. Reads left by two trimming routes stop
`evaluate`.

The genome-wide synthetic set is drawn at the filter's minimum read length (#500): the declaration's
`config.filter.trim.cutadapt.minimum_length`, otherwise its kit's in `libprep.config` (the set runs the cutadapt route).
No read of the truth, of the reads another miRNA's loci could give, of the decoys or of the random reads is shorter. `prepare` records the
minimum as `synthetic.min_length` in the project's `benchmark/dataset.yaml`; a declaration cannot set it. That record
is the length the truth was drawn at, not the one a later run filters at. A project drawn at 15 nt and run with the
filter at 16 keeps 15-nt reads in its truth that no route sees (3′ ends shortened and 5′ ends moved on 16- and 17-nt
miRNAs); they count as not recovered for every route alike. A project drawn at 16 has none. Projects prepared before
#500 have no `min_length` and were drawn at 15; `summary` reads them so and groups their seeds with draws recorded at 15.
At the same minimum and seed the draw is the same as before #500. `prepare` of this set reads the kit with configmaker,
as the workflow does, so it runs in an environment that has it (the workflow's own, or the repository's `.dev/ci`).

**Runs from SRA instead of ENA.** ENA serves a client one slow connection and may cut or refuse transfers (`fill`
resumes a cut one). NCBI's copy of the Sequence Read Archive on AWS takes many requests at once:

```bash
python smallrna/benchmark/benchmark.py fill --dataset <declaration> --ext-dir "$EXT_DIR" --source sra \
    --fasterq-dump "apptainer exec --contain --bind $EXT_DIR docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1 fasterq-dump"
```

- **What a run needs:** `sra: {md5: ..., bytes: ...}` in its entry, the SRA file's checksum and size from NCBI's
  locator (`https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=<run>`: the file of type `sra` named as the run).
  A run without it stops the fill; runs already in place are skipped, whatever their source.
- **What happens:** the SRA file is downloaded in parts (`--connections` at a time, `--jobs` runs at a time),
  checked against its MD5, converted with `fasterq-dump` (`--fasterq-dump`: the command; sra-tools 3.2.1, as 3.4.1's
  image crashes) and written as `fastq/<run>.fastq.gz`, with `fastq/<run>.sra.json` recording the source, the
  checksums and the tool version.
- **How it is checked:** the reads are those of ENA's file in the same order, with other read names, so the file's
  bytes and MD5 differ from the declared ENA `md5`. The check is the number of reads (`reads`) and, when the run
  declares it, `content_md5`: the MD5 over the sequence and quality lines, which is the same for a file from either
  archive (`fill.file_content` computes it for a file in place). A fill without a declared `content_md5` prints the
  value to declare.
- **Limits:** single-end runs only. `prepare` cannot check a declared lane against read names from SRA and says so.

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
| `gse317995_nextflex_v4` | human plasma, NEXTflex v4 (GEO GSE317995): all 30 samples, 10 in each of three groups (#448); `prepare --per-group 2` takes the two lowest-numbered of each, the quick run | 4 per sample: the 4 lanes of one NextSeq 500 run |
| `gse77989_hct116_knockouts` | HCT116 cells, TruSeq Small RNA (GEO GSE77989, Kim et al. 2016): parental cells (3), DROSHA (2), XPO5 (2) and DICER (2) knockouts, in two sets of libraries | 1 per sample (MiSeq) |
| `gse271895_nextflex_v4_lps` | human fibroblasts, NEXTflex v4 (GEO GSE271895, McLaurin et al. 2025): a scenario of `gse271895_nextflex_v4` grouped by drug, LPS (11) against DMSO (12) with the siRNA arm as a covariate and per arm (`within: sirna`), and null splits; the 55,677-read library is left out (#414) | 1 per sample (NextSeq 2000) |

A biological dataset with knockout cells declares `knockout:` (#452): the `control` group, the knockout whose
surviving miRNAs are listed (`survival`), and optionally `within`, a field of the samples that pairs each knockout sample
with the controls of its own level. GSE77989 has `within: set`: the DROSHA knockout was sequenced with two parental
libraries (set 1), the XPO5 and DICER knockouts with one (set 2). Its knockout clones are 2 libraries each, shallow, from
one cell line after gel selection; XPO5 knockout lowers miRNAs only partly, so it is reported but not used for the list.

`knockout.expected` (#478) names a list of the miRNAs made without DROSHA, a file next to the declaration that `prepare`
copies into the project's `benchmark/`. GSE77989's `knockout_expected_hsa.tsv` is written by `knockout_expected.py` from
two sources, each recorded in its header:
- **Kim et al. 2016, Fig. 3B:** the study's own list of the miRNAs its DROSHA knockout keeps, among the 200 most
  abundant: the 5' capped miR-320a-3p, miR-320b-3p, miR-484 (miRBase's hsa-miR-484-5p), miR-7706-3p and miR-3615-3p,
  and the mirtron miR-877-5p. Its miR-1254-5p, from an Alu element, is not in miRBase 23.1.
- **mirtronDB version 2:** every miRBase mature that lies exactly inside one of its 607 human mirtron precursors (469
  matures), with the PubMed IDs of the papers behind them. Most are predictions from small RNA data, not tested in
  these cells.

With `other_rna` declared, the flagged miRNAs (inside another RNA class) are expected too.

A biological dataset may declare `other_rna:` (#415), the sources of molecules of other RNA classes: `refseq_rrna`,
`gtrnadb` and `ensembl_ncrna`, each with url and sha256, the same files as the synthetic decoys, filled into
`<data>/other_rna/`. Its miRBase entries whose mature lies exactly inside such a molecule (rRNA, tRNA, Y RNA, snoRNA,
other ncRNA) are flagged (`other_rna.tsv`); a fragment that overlaps a mature only partly is not. The report and
metrics.tsv then give, per method:
- `flagged_reads`: the median share of its miRNA reads on flagged entries;
- `flagged_detected` and `detected_unflagged`: the median detected miRNAs among the flagged and without them, so
  `detected` can be compared without the inflation (#385); only `detected_unflagged` has an interval;
- `read_in_other_class`: the share of its miRNA reads (isomiR output, #454) whose read sequence lies exactly inside
  another class's molecule.

Both the flag and `read_in_other_class` are lower bounds: they cover five classes (no mRNA, introns or repeats) and
exact matches only, so a read with a non-templated end or a sequencing error does not match.

With knockouts, `knockout.tsv` names each tested miRNA's other class, and `knockout_surviving_flagged` counts the
surviving miRNAs that are flagged. Both human biological declarations have `other_rna`.

GSE317995 is the only public NEXTflex v4 dataset found. Its lanes are one library each, so they replicate sequencing
only; library-preparation reproducibility needs separate libraries of one sample (e.g. the Giraldez 2018 plasma pool
made in several labs). `blockers: [blood miRNA]` is inferred, not recorded: miR-486-5p and miR-451a are nearly absent.

**Tissue mixtures** (`kind: mixture`, #451): biological runs of a known structure. Total RNA of two tissues
(`mixture: {components: [first, second]}`) of each subject is mixed at known fractions. A declaration lists `samples`
as a biological one does, each with:
- `subject`: the animal or donor whose two tissues were mixed;
- `fraction`: the share of the second tissue in the total RNA, from 0 to 1; 0 and 1 are the subject's pure samples;
- `replicate` (optional, from 2): a further library of the same subject and fraction, a technical replicate.

A subject's series is its two pure samples and its mixtures. A sample of a subject without both pure samples, such as
a pool of several subjects' RNA, is outside the series and not scored. The composition of no sample is known, but each
mixture is the weighted sum of its subject's pure samples, so every method is scored against its own counts (Metrics,
below). The project is prepared and quantified as a biological one, without `--lanes` and without `--notemp`: the
evaluation reads only the count tables.

| Declaration | Samples | Runs |
|---|---|---|
| `srp697344_liver_spleen_nebnext` | mouse liver and spleen of five mice, each at 0, 20, 40, 60, 80 and 100% spleen, with a technical replicate at each fraction that has spleen and a pool of the five livers; NEBNext (SRA SRP697344; Aparicio-Puerta et al., bioRxiv 2026, doi:10.64898/2026.05.08.723923, a preprint) | 36, one per library, 17–40M reads |

The study's own truth table is not used: it is derived from one tool's counts and would favour methods that count
like it. The authors note pancreas RNA in mouse 3; the report gives the error per subject, so an outlying subject
shows.

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
  reads in the project (as `miR.Counts.csv`); a read under `exact miRNA` is labelled canonical. miRge3's exact class
  holds reads matching its library's sequences, not only the miRBase mature: 82% are the mature sequence on the
  genome-wide synthetic set, the rest mostly 3′ shorter or longer, so its label metrics describe that class.
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
- **Differential expression** (#447; `de: {}` and `evaluate --de`): the tests of the ratiometric pools (Pools below)
  on the mixed samples, so the two datasets can be compared. The group the designed fold changes raise is tested
  against the other, with the workflow's library sizes (all the miRNA reads a method reported; the truth: its miRNA
  reads). Each row tests the miRNAs with at least 10 reads on average in its own counts, as an analysis of them would,
  so reads spilled onto a relative of a changing miRNA can make a false call. The scored miRNAs are those designed to
  change with at least 10 truth reads on average, and every other miRNA of the truth or reported by a method that is
  designed not to change. The metrics are the pools', with the same prefixes: `called`, `false_calls`, `fdr`,
  `fc_bias` and `offset` against the designed changes, and against the truth's fits `own_false_calls`,
  `discordant_vs_truth`, `fc_shift_vs_truth` and `fc_error_vs_truth`. A second set adds the nuisance factor's true
  sample scores as a covariate of the GLM test and voom; its metrics end in `.adjusted` (the exact test takes no
  covariates, so its values there equal the unadjusted ones). The report gives the share called per designed fold
  change, and the share of each method's miRNAs whose fractional counts were rounded for the tests.
- **Differential expression, the built-in check** (`de_*`, without a container): median-of-ratios size factors,
  log2(normalised + 1), a linear model per miRNA with the group (and, `_adjusted`, the true nuisance scores), a
  moderated t (the residual variance shrunk towards the median, 4 prior degrees of freedom), Benjamini-Hochberg at 5%,
  on the miRNAs with at least 10 reads on average in a row's own counts (`de_tested`). Reported as `de_sensitivity`,
  `de_fdr` (calls on any miRNA not designed DE, absent ones included; the report lists them) and `de_fc_error` (median
  |observed - true log2 fold change| of the true DE miRNAs). It differs from the tests above in normalisation, model
  and scoring, so compare methods within it, not its values with the pools'.
- **Structure:** `structure_group_r2` and `structure_nuisance_r2` (each true factor regressed on the top two principal
  components) and `coexpression` (agreement of the miRNA-miRNA correlations).
- **Probe samples:** `probe_<class>`, the recovery of each read class, control subtracted. A generated read can still be
  a variant of another mature beyond these classes (a 2-nt 3′ trim of miR-320d-3p is the 5′ 17 nt of the miR-320
  family). The truth keeps its origin, which no counting method can know, so a method may count such reads on a
  relative.
- **Decoys** (#412 PR 2): `decoy_as_mirna_<class>-<source>`, the share of a decoy probe's reads that a method counts
  on miRNAs other than the carrier, control subtracted; `decoy_as_mirna` pools the units. Every such read is false.
  `probe_embedded_<class>` is the recovery of the miRBase matures inside that class's molecules, so 1 − recovery is the
  share lost to the other class. `overlap_on_mirna_<class>` and `overlap_as_mirna_<class>`: the share of the overlap
  fragments counted on the mature they overlap and on any miRNA (not false positives, as miRBase has the mature). The
  report lists the miRNAs decoys land on. Decoys in the mixed samples act on
  `precision` and `misassigned`. There are no intervals, as for the probes.
- **isomiRs** (#454; mixed samples; intervals resample the miRNAs, #459): a true read counts as reported when the method reports its own
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

  **Normalisation arms** (#453): the whole design is tested again in two other normalisations, for the pools and the
  genome-wide synthetic set, and the report's Normalisation arms section shows them beside the workflow's:
  - `fixed` (folder `de/arm/fixed`, metrics `.fixed`): every row's library size is the sample's filtered reads, the
    cutadapt route's kept reads, the same for every method, so what a method counts as miRNA no longer changes it;
  - `tmm` (folder `de/arm/tmm`, metrics `.tmm`): edgeR's TMM factors (`calcNormFactors`) on the workflow's library
    sizes, for all three tests (voom on the normalised sizes).
  `false_calls_against` (all DE sets, #453): the share of the unchanged features called in the direction opposite to
  the known changes, when at least two thirds of them go one way (empty otherwise, as in the pools and the balanced
  synthetic set).
  The main metrics stay the workflow's. An `offset` that differs between the arms is the normalisation's, not the
  method's. Biological datasets have no arms: no change is known to score them against. The `fixed` arm reads the
  cutadapt route's reports (the kits' default), so `evaluate` of a pool or synthetic project that declares `de:` on
  another filter route stops and names the missing reports. On the pools `false_calls` is a share of their 25 1:1
  sequences and, with dispersions near 0.002, follows the offset, so read the arms by their `offset`.

  With `de: {within: lab}` (#444) each lab's runs are also tested on their own with `~ group`, as replicate benchmarks
  of a single-lab experiment; each lab needs at least two runs per pool. Their metrics carry the lab, e.g.
  `edger_fc_error_vs_oracle.Lab1`. They are in `metrics.tsv`, `metrics_ci.tsv` and `differences.tsv`, whose paired
  differences show whether a difference between methods holds in each lab. The report shows their values in a section
  of their own, outside the main table and the uncertainty section.
- **isomiRs** (#454; intervals resample the scored sequences and the runs, #459). Against the oracle, per run: `exact_sequence_recovered` and
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

- **Differential expression** (#447; `de:` and `evaluate --de`): the pools' tests between two groups of samples, each
  method on the miRNAs with at least 10 reads on average in its own counts, identical sequences as one. The report
  gives what each method tests and calls up and down per test. No change is known, so nothing is scored and
  `metrics.tsv` gets no DE metrics.

- **Knockouts** (#452; `knockout:`): canonical miRNAs need DROSHA and DICER, so their reads should nearly vanish in
  those knockouts. Reads per million filtered reads (the cutadapt route's kept reads, not miRNA reads, which the
  knockout changes); each knockout sample against the control samples of its `within` level, merged samples only.
  - `knockout_mirna_share`: the survival knockout's miRNA reads per million filtered reads over its controls' (mean
    over its samples); the report has it for every knockout group, and the controls' miRNA reads per million.
  - The ratios are upper bounds: they are shares of the filtered reads, so when the miRNAs go, the other RNAs take a
    larger share. In GSE77989's DROSHA knockout (miRNAs 60-71% of the parental filtered reads, 3-7% in the knockout)
    an RNA whose amount did not change rises about 3x, so a miRNA kept at 25% has fallen to about 8%; the report
    gives the factor per dataset. The routes share the denominator, so they compare.
  - Per miRNA with at least 10 reads on average in the controls: the knockout's reads per million over the controls';
    it survives above 25%. `knockout_surviving` counts those of the survival knockout; the report lists them with
    each method's ratio, those every method lists, and per method those only it lists (`only`) with their reads.
    Surviving miRNAs common to every method are biology (e.g. mirtrons) or annotation; those that differ are
    quantification. `knockout.tsv` has every tested miRNA of every knockout group.
  - With `knockout.expected` (#478): the miRNAs expected to survive, made without DROSHA (the list) or inside another RNA
    class (`other_rna`). `knockout.tsv` gives each one's reason.
    - The report counts per method the expected miRNAs it tests and those that survive (`knockout_surviving_expected`),
      and names those tested but not surviving.
    - It splits `only` into the expected, those another method tests without listing them (it measures them lower),
      and those no other method tests. That split depends on which methods run; `knockout_surviving_unexpected`, the
      surviving miRNAs not expected, does not.
    - It names the miRNAs surviving in every method (two or more) that are not expected: candidates to check for the
      list rather than one method's false counts. `expected` in `knockout.tsv` is filled for the survival knockout only.
- **Accessions.** miRNAs with the same mature sequence count as one (the lowest accession), since the methods credit
  their reads to different members (miRge3 to its merge group's first name; UNITAS and mirflowz to either, e.g.
  hsa-miR-199a-3p and hsa-miR-199b-3p). Reads on names without a miRBase accession are left out of every metric and
  reported as `unplaced`, with the names that have most.
- **Per method:** miRNA reads, `unplaced`, miRNAs detected per sample and in every sample, and the share of the 10 most abundant.
- **isomiRs** (#454), a description: `isomir_share` (the share of a method's miRNA reads that are not the mature
  sequence; its interval resamples the samples, #459) and `isomirs_detected` (sequences with at least 10 reads over all
  samples; no interval, as `detected_all_samples`); the report adds each method's reads by kind and the overlap of
  detected isomiRs between methods (Jaccard).
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

**Tissue mixtures** (`mixture.py`, #451), per method, against its own counts of the pure samples:

- **The model.** With `fraction` f of the total RNA from the second tissue, the share of a mixture's miRNA reads from
  it is w = f k / (f k + 1 − f). The content ratio k is the second tissue's miRNA reads per unit of RNA over the
  first's, as the method counts them; it is fitted per method and subject, since each subject's tissue samples have
  their own RNA content (each mixture's weight by least squares of its proportions on the subject's two pure
  profiles, turned into k; the median over the subject's mixtures). A per-miRNA bias that is the same in every sample
  changes k and nothing else, so the kit's bias does not enter the scores.
- **`linearity_error`**: the median |log2(observed / predicted)| over the miRNA-sample pairs of the mixtures with at
  least 10 predicted reads; **`within_1_5fold`**: the share within 1.5-fold. The report also gives the error per
  subject and by predicted reads.
- **`monotonic`**: of the tissue-specific miRNAs (at least 4-fold apart between a subject's pure samples, 10 predicted
  reads in every mixture), the share of subject series that never move against the trend over the fractions.
- **`replicate_error`**: the median |log2 ratio| of the proportions of technical replicates, miRNAs with at least 10
  reads on average in the pair. Two libraries of one RNA differ by this much, so it is the floor for
  `linearity_error`.
- **`content_ratio`** (the median of the subjects' k) and **`content_ratio_spread`** (the median |log2| deviation of
  the subjects from it): a description of the fit, not a score. The report gives each subject's k.
- **What it shows.** Linear and consistent counting across samples. It does not show correct assignment: a method that
  puts a read on the wrong miRNA in every sample still passes. The pools and the synthetic dataset score assignment.
- **Differential expression** (`de: {}` and `evaluate --de`), with the pools' tests. The GLM test and voom pair the
  samples by subject (`~ subject + group`); the exact test takes no covariates. Library sizes are a method's reads on
  miRNAs with an accession, the total its proportions and so its expected changes use:
  - the pure samples, second tissue against first: what each method tests and calls, as for biological samples; not
    scored, since the pure samples define the expectation;
  - each pair of the mixtures' fractions: a method's expected change of a miRNA is the mean over the subjects of
    log2(prediction at one fraction / prediction at the other), from pure samples that are in none of these
    contrasts. Over the miRNAs a method tests with an expected change of at least 0.2 log2, all contrasts together,
    prefixed `edger_`, `edger_glm_` and `voom_`: `called` (expected at least 2-fold apart and called in the expected
    direction), `wrong_direction` (the share of its calls against it), `fc_bias` (as for the pools),
    `called_shared` (`called` over the miRNAs every method expects at least 2-fold apart in the same direction, so
    the methods are compared on the same miRNAs) and `unchanged_called` (the share called of the miRNAs with a
    smaller expected change than 0.2 log2; hardly any miRNA is truly unchanged between two mixtures, so it compares
    the methods and is not a false-call rate).
- **Intervals** resample the miRNAs and the subjects on the fitted content ratios (`replicate_error`: the miRNAs
  only); the DE metrics have none.

On the mouse chr19 dataset, UNITAS and miRge3 count every miRNA exactly. mirflowz counts a third of
`mmu-miR-146b-5p`, whose sequence occurs three times on chr19: it divides a read that maps equally well to several loci
between them, while UNITAS and miRge3 match by sequence.
