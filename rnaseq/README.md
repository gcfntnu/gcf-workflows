# rnaseq

**Status: active.** Maintained and covered by a CI dry run. Maintenance status does not mean every kit, route or output listed here has been scientifically validated.

Bulk RNA sequencing: read QC and trimming, genome alignment, transcript and gene quantification, QC metrics and exploratory figures.
Entry point: [`rnaseq.smk`](rnaseq.smk). For how BFQ runs a workflow and what it delivers, see the [README](../README.md#how-a-workflow-runs-under-bfq).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: rnaseq`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts |
| --- | --- |
| Lexogen SENSE mRNA-Seq Library Prep Kit V2 | SE, PE |
| Lexogen SENSE Total RNA-Seq Library Prep Kit (w/RiboCop rRNA Depletion Kit V1.2) | SE, PE |
| CORALL Total RNA-Seq Library Prep Kit (w/RiboCop rRNA Depletion Kit V1.2) | SE, PE |
| Lexogen CORALL mRNA-Seq Library Prep Kit | SE, PE |
| Illumina TruSeq Stranded mRNA Library Prep | SE, PE |
| Illumina Stranded mRNA Prep | SE, PE |
| Illumina Stranded Total mRNA Prep with Ribo-Zero Plus | SE, PE |
| Illumina TruSeq Stranded Total RNA Library Prep (Globin) | SE, PE |

The kit entries set the strand orientation, adapters, trimming parameters, read-length adjustment and `db.reference_db: ensembl`.

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| filter | `rules/filter.smk`, shared `common/` rules | FASTQ merging and fastp trimming |
| assembly | `rules/assembly.smk` | Transcript and genome assembly rules, loaded only when an assembler is not `skip` (the default) |
| align | `rules/align/star.smk`, `hisat2.smk` | STAR by default (`align.genome.aligner`); HISAT2 available |
| quant | `rules/quant/salmon.smk`, `featurecounts.smk`, `tximport.smk` | Salmon by default (`quant.method`), featureCounts, tximport |
| qc | `rules/qc/fastq.smk`, `bam.smk`, `quant.smk`; shared `common/` rules | FastQC and Kraken2 screening (shared); BAM (Picard RNA metrics, RSeQC TIN) and quantification QC |
| bfq | `rules/bfq/` (levels 1–3) | Collect delivery products into `bfq/` |

**Inactive source.** Files present in the source tree are not necessarily loaded. The entry point does not load these, so configuration cannot select them:

- `rules/analysis/` (`preprocess`, `diff_expr`, `co-exprs`): `rules/analysis.smk` includes nothing.
- `rules/quant/` alternative splicing, CNVkit, HLA typing, Rockhopper, T-cell repertoire and variant calling.

## Outputs in `bfq/`

`data/tmp/rnaseq/bfq/`, as planned by the fixture's dry run (single-end, Salmon, STAR). Paired-end data has `_R1`/`_R2` FastQC files. Other methods change the set.

```text
multiqc_<project>.html, .multiqc_config.yaml
align/<sample>.sorted.bam, .sorted.bai, .sorted.bam.md5
exprs/gene_counts.tsv, gene_tpm.tsv, gene_vst.tsv, gene_info.tsv
exprs/transcript_counts.tsv, transcript_tpm.tsv, transcript_vst.tsv, transcript_info.tsv, tx_salmon.rds
figs/gene_biotypes_mqc.yaml, gene_high_mqc.yaml, pca_mqc.yaml
logs/krona_all_samples_kraken.html
logs/<sample>/<sample>.fastp.html, .fastp.json, _fastqc.html, _fastqc.zip
logs/<sample>/<sample>.kraken.kreport, .kraken.out
logs/<sample>/<sample>.Log.final.out, .rnaseq.metrics, .tin.summary.txt
logs/<sample>/aux_info/meta_info.json
```

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | `salmon` | [rnaseq.config](rnaseq.config) |
| `align.genome.aligner` | `star` | [rnaseq.config](rnaseq.config) |
| `align.star.args`, `align.star.guided_index` | STAR options | `rnaseq.config` |
| `db.reference_db`, `db.ensembl.release` | `ensembl`, release in `main.config` | Kit entry, `main.config` |
| `read_orientation`, `adapter`, `filter.trim` | Per kit | Kit entry |
| `assembly.*` | `skip`; `merge.type: fastq` | `rnaseq.config` (preserves the previously effective `fastq` default) |
| `qc.bam_subset` | `housekeeping_genes` | `rnaseq.config` |

Precedence between these sources is described in [docs/configuration.md](../docs/configuration.md#effective-precedence).

## Tests

- `rnaseq/.tests`: synthetic dry-run fixture, single-end, kit `Lexogen SENSE Total RNA-Seq …`, `homo_sapiens`. CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-rna-seq.yaml`).
- The [standalone recipe](../README.md#running-standalone) dry-runs rnaseq with a paired-end synthetic project.

A dry run shows only that the job graph can be built. See [what the checks prove](../docs/development.md#what-the-checks-prove).

## More documentation

- [`methods.txt`](methods.txt): methods text for reports. It is not read by the workflow.
- [Architecture](../docs/architecture.md), [configuration](../docs/configuration.md), [reference data](../gcfdb/README.md)

## Credits

Written by Arnar Flatberg ([flatberg](https://github.com/flatberg)) at the Genomics Core Facility, [Norwegian University of Science and Technology, NTNU](http://www.ntnu.no). Licence: [LICENSE](LICENSE).
