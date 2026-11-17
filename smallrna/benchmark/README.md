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
3. **Run** `snakemake quant_all` in the project, as for any project on that machine (containers, cores, resources).
4. **Evaluate.** `benchmark.py evaluate --project <project>` writes `<project>/benchmark/`: `common.tsv` (every
   method's counts in one table), `per_mirna.tsv`, `metrics.tsv` and `report.md`.

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
`EXT_DIR/<data>/` (written under a temporary name, renamed after the checksum matches); it is the only step that
writes to `EXT_DIR`:

```bash
python smallrna/benchmark/benchmark.py fill --dataset smallrna/benchmark/datasets/giraldez2018_equimolar_truseq.yaml \
    --ext-dir "$EXT_DIR"
```

`prepare` links the runs into the project, or copies the first `--reads N` reads of each for a quick run; `--methods`
picks a subset of the declared methods (mirflowz on human needs the full human reference). Run `quant_all` with
`--notemp`: the oracle count reads the workflow's filtered reads, those of the project's filter route (`filter/cutadapt/`
with the kit default; the spike-in prefilter's with spike-ins). Reads left by two routes stop `evaluate`.

## Methods

`adapters.py` turns each method's native output into the common table: method, sample, the method's miRNA name, the
miRBase accession and the count (canonical reads and isomiRs together). Matching is by accession, because the methods
name miRNAs differently: UNITAS drops the species prefix, miRge3 (the workflow's own `miR.Counts.csv`, built from
miRge3's `mapped.csv` with the project-level exact-read rule, #349) reports merge groups (`name1/name2`, counted on the
first member), mirflowz names miRNAs after their locus (`mmu-miR-194-2-5p`) and reports isomiRs as separate rows. A new
method needs one function in `ADAPTERS`, next to its smallrna route.

## Metrics

**Synthetic data.** For truth as reads per miRNA, per method: `recovered` (share of the true reads counted on the right
miRNA), `exact` (share of miRNAs counted exactly), `misassigned` (reads on miRNAs not in the truth) and `unplaced`
(reads on names without a miRBase accession).

**Pools** (`pool.py`):

- **Truth groups.** As in the study, sequences of 16-25 nt with only a 5' phosphate and unique in the pool are
  analysed, grouped against the organism's miRBase mature sequences: `exact` (the sequence of one or more mature
  miRNAs), `shifted` (a mature miRNA with up to 3 nt more or less at either end: expected as an isomiR), `foreign`
  (other species, viral, revised or non-miRNA sequences) and `excluded`. Against human miRBase 23.1 the equimolar pool
  has 407 `exact` sequences, the ratiometric pools 211.
- **Oracle.** The workflow's filtered reads counted by exact match to each pool sequence: what a sample contains after
  the kit's bias. Methods are scored against it on the `exact` group, so a kit's per-sequence bias (up to 100-fold) does
  not hide the differences between methods.
- **Equimolar pool**, per method and for the oracle: `detected`, `cv` and `qcd` (quartile coefficient of dispersion)
  across the group; per method against the oracle: `median_log2_vs_oracle` (of (method + 1) / (oracle + 1)),
  `within_2fold`, `beyond_10fold`, `spearman_vs_oracle`, and `outside` (share of the method's miRNA reads on miRNAs
  not in the pool, i.e. assigned from foreign sequences).
- **Ratiometric pools**: counts per million of the group, averaged per pool, give an observed log2(A/B) per sequence,
  centred on the sequences designed 1:1. Against the design: `median_abs_error`, `spearman_vs_design`, `called`
  (designed at least 2-fold apart and observed so, in the right direction) and `false_calls` (designed less than 2-fold
  apart, observed at least 2-fold apart). A fold change cancels the per-sequence kit bias, so what differs between
  methods is quantification. The report also lists the median observed log2(A/B) per designed ratio.

On the mouse chr19 dataset, UNITAS and miRge3 count every miRNA exactly. mirflowz counts a third of
`mmu-miR-146b-5p`, whose sequence occurs three times on chr19: it divides a read that maps equally well to several loci
between them, while UNITAS and miRge3 match by sequence.
