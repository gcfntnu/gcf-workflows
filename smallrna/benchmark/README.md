# Small-RNA quantification benchmark

Compares the smallrna quantification methods (`quant.method`: UNITAS, miRge3, mirflowz) on datasets with a known
truth, and makes adding a method or a dataset a small change ([#338](https://github.com/gcfntnu/gcf-workflows/issues/338)).
It runs the workflow's own `quant_all` routes, never the tools directly, so it measures what a project would get.
It is not run in CI; its own code is tested in `smallrna/.tests/test_benchmark.py`.

## Steps

1. **Reference.** A reduced reference (a few chromosomes) lets mirflowz run on a workstation. It goes to
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

Declarations live in `datasets/`. A synthetic declaration lists miRNAs with reads per sample (`mirnas`), miRNAs drawn
from those annotated on the reference (`random_mirnas`), isomiRs (`isomirs`: `3p_trim_1`, `3p_add_A`, `5p_trim_1`) and
random reads; reads are generated from a fixed seed with the kit's adapter (and `umi` random nucleotides at each end
for kits that use them), so the workflow's own filtering runs as usual. Public datasets with a known composition will
be added with a fill step that downloads them into `EXT_DIR/benchmark/smallrna/<dataset>/` once, checked against
recorded checksums; benchmark runs only read them.

## Methods

`adapters.py` turns each method's native output into the common table: method, sample, the method's miRNA name, the
miRBase accession and the count (canonical reads and isomiRs together). Matching is by accession, because the methods
name miRNAs differently: UNITAS drops the species prefix, miRge3 reports merge groups (`name1/name2`, counted on the
first member), mirflowz names miRNAs after their locus (`mmu-miR-194-2-5p`) and reports isomiRs as separate rows. A new
method needs one function in `ADAPTERS`, next to its smallrna route.

## Metrics

For truth as reads per miRNA (synthetic data), per method: `recovered` (share of the true reads counted on the right
miRNA), `exact` (share of miRNAs counted exactly), `misassigned` (reads on miRNAs not in the truth) and `unplaced`
(reads on names without a miRBase accession). For an equimolar pool: share of pool miRNAs detected, coefficient of
variation, and share within two-fold of the median.

On the mouse chr19 dataset, UNITAS and miRge3 count every miRNA exactly. mirflowz counts a third of
`mmu-miR-146b-5p`, whose sequence occurs three times on chr19: it divides a read that maps equally well to several loci
between them, while UNITAS and miRge3 match by sequence.
