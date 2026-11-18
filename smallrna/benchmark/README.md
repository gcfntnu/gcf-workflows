# Small-RNA quantification benchmark

Compares the smallrna quantification methods (`quant.method`: UNITAS, miRge3, mirflowz) on datasets with a known
truth, and makes adding a method or a dataset a small change ([#338](https://github.com/gcfntnu/gcf-workflows/issues/338)).
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
   method's counts in one table), `per_mirna.tsv` (pools: `per_sequence.tsv` and `filter_stages.tsv`), `metrics.tsv`
   and `report.md`; with UNITAS, `unitas_names.tsv` (how each UNITAS name was placed, below). UNITAS names are placed
   with the project's UNITAS reference dump (`db.unitas` under `EXT_DIR`, read only; `--unitas-ref` names another).

## Datasets

Declarations live in `datasets/`.

**Synthetic** (`kind: synthetic`): miRNAs with reads per sample (`mirnas`), miRNAs drawn from those annotated on the
reference (`random_mirnas`), isomiRs (`isomirs`: `3p_trim_1`, `3p_add_A`, `5p_trim_1`) and random reads; reads are
generated from a fixed seed with the kit's adapter (and `umi` random nucleotides at each end for kits that use them),
so the workflow's own filtering runs as usual. Not stored.

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
picks a subset of the declared methods (mirflowz on human needs the full human reference). Run `quant_all` with
`--notemp`: the oracle count reads the workflow's filtered reads, those of the project's filter route (`filter/cutadapt/`
with the kit default; the spike-in prefilter's when the project's `filter.spikein.ref` lists sets removed by Bowtie,
although the trimmed reads it starts from remain; the trimmed reads with native spike-ins, which UNITAS and miRge3
read), and the filter-route count reads those after each of the route's passes. Reads left by two trimming routes stop
`evaluate`.

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
TEsmall reports mature miRNAs by name (precursor counts left out), and tiny-count reports one row per locus by
miRBase feature ID (`MIMAT0000416_1`), summed per accession. A new
method needs one function in `ADAPTERS`, next to its smallrna route.

A UNITAS name is placed by name when it is a mature miRNA of the organism, else by hairpin and arm: the hairpin of that
name in the UNITAS dump (`miRBase_hairpin_<Species>.fas`), and the organism's mature miRNA on that arm, the hairpin
family's (`miR-16-1-5p`: hsa-miR-16-5p) or the only one. Other names stay unplaced: no mature miRNA on that arm
(human miRBase has only hsa-miR-325-5p, so `miR-325-3p`), several of other families, or no such hairpin. Stripping
the locus number is not enough: hsa-miR-181a-2-3p is its own mature miRNA. On the Giraldez equimolar pool (TruSeq) this
places 4% of UNITAS's miRNA reads that had no accession and leaves 0.2% unplaced; the reads placed are those the other
methods count on the same miRNAs (often a pool sequence one or two bases shorter).

## Metrics

**Synthetic data.** For truth as reads per miRNA, per method: `recovered` (share of the true reads counted on the right
miRNA), `exact` (share of miRNAs counted exactly, to within 1e-6 for the fractional counts of divided multimappers), `misassigned` (reads on miRNAs not in the truth) and `unplaced`
(reads on names without a miRBase accession).

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
- **Between methods**, per pair (medians over samples): Spearman of log2(CPM + 1) over the miRNAs either detects
  and over those both detect, median |log2 CPM ratio| over those both detect, overlap of the 100 most abundant, and miRNAs only one detects
  (`only_a` and `only_b` for the first and second method of the pair).
- **Blocked miRNAs:** CPM of miR-451a, miR-486-5p and miR-92a-3p, the targets of the NEXTflex blood miRNA blockers.
- **Filtering:** reads through the cutadapt route's three passes (3' adapter found, 5' adapter discarded, too short,
  N, expected errors, artefacts discarded, kept), from their JSON reports, which must chain (each pass read what the
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
