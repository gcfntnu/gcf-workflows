# microbiome

> **Unmaintained.** This workflow is kept in the repository but is not maintained. It may not run, and its results are not validated. Its CI dry run still passes, which shows only that the job graph can be built.

Amplicon (16S/ITS) microbiome analysis based on QIIME2: primer and region demultiplexing, denoising, taxonomic classification, phylogeny, and phyloseq-based analysis notebooks.
Entry point: [`microbiome.smk`](microbiome.smk).

## Library-prep kits

Kits whose [libprep.config](../libprep.config) entry sets `workflow: microbiome`. `libprep.config` is authoritative; this list reflects it at the time of writing.

| Kit | Layouts | `db.reference_db` |
| --- | --- | --- |
| 16S Metagenomic Sequencing Library Prep | PE | `silva` |
| ITS Low Input GCF Custom | PE | `unite` |
| ITS2 Low Input GCF Custom | PE | `unite` |

The `QIAseq 16S ITS Region Panels` kit is routed to [metagenome](../metagenome/README.md). The microbiome test fixture uses that kit with `workflow: microbiome` set explicitly.

## Stages

| Stage | Rule files | Notes |
| --- | --- | --- |
| reference_db | `rules/reference_db.smk` | SILVA, Greengenes or UNITE, by `db.reference_db` |
| filter | `rules/filter/`, shared `common/` rules | Cutadapt primer trimming and QIAseq region demultiplexing; fastp |
| qc | `rules/qc/fastq.smk`, shared `common/` rules | FastQC and FastQ Screen |
| quant | `rules/quant/qiime2_py.smk`, `qiime2_filter.smk`, `qiime2_diversity.smk`, `qiime2_picrust.smk` | QIIME2 (`quant.method: qiime2`): denoising, filtering, diversity, PICRUSt |
| analysis | `rules/analysis/phyloseq.smk`, `picrust.smk` | phyloseq notebooks, PICRUSt |
| bfq | `rules/bfq/` (levels 1–2) | Collect delivery products into `bfq/` |

**Inactive source.** Files present in the source tree are not necessarily loaded. The entry point does not load `rules/qc/bam.smk`, `rules/qc/quant.smk`, `rules/quant/qiime2.smk`, `rules/analysis/qiime2_gneiss.smk` or `qiime2_longitudal.smk`.

## Known problems

- The `phyloseq_filter` notebook sources `microfiltR` code that is not in the repository ([#231](https://github.com/gcfntnu/gcf-workflows/issues/231)).
- QIIME2 report support is deprecated for the MultiQC migration ([gcf-multiqc](https://github.com/gcfntnu/gcf-multiqc)).
- `ref_dir`/`GCF_REFDIR` resolve in a different order here than in most workflows ([#240](https://github.com/gcfntnu/gcf-workflows/issues/240)).

## Outputs in `bfq/`

`data/tmp/microbiome/bfq/`, as planned by the fixture's dry run (QIAseq regions, paired-end):

```text
multiqc_<project>.html, .multiqc_config.yaml, aggr.done
exprs/table.biom, table.tsv, dna-sequences.fasta, feature_info.tsv, samples_info.tsv, tree.nwk, physeq.rds
figs/qiaseq_regions_mqc.yaml
logs/deicode/ordination.txt, ordination_samples.txt
logs/taxa/level-1.csv … level-7.csv
logs/<sample>/<sample>.fastp.html, .fastp.json, _R1_fastqc.*, _R2_fastqc.*, _screen.txt
qiime2/
region_demultiplexed_fastq/<sample>_R1.fastq.gz, <sample>_R2.fastq.gz
```

## Key configuration

| Key | Default | Set in |
| --- | --- | --- |
| `quant.method` | `qiime2` | [microbiome.config](microbiome.config), kit entry |
| `db.reference_db`, `db.primers`, `db.classifier_level` | Per kit | Kit entry |
| `db.silva.*`, `db.unite.*`, `db.greengenes.*` | Releases | [main.config](../main.config) |

## Tests

`microbiome/.tests`: synthetic dry-run fixture (four samples, two groups). CI runs `snakemake -n -p --cores 1 multiqc_report` there (`test-microbiome.yaml`).

## More documentation

- [docs/qiaseq-demultiplexing.md](../docs/qiaseq-demultiplexing.md): QIAseq region demultiplexing, shared with metagenome
- [Architecture](../docs/architecture.md), [configuration](../docs/configuration.md)

## Credits

Written by Arnar Flatberg ([flatberg](https://github.com/flatberg)) at the Genomics Core Facility, [Norwegian University of Science and Technology, NTNU](http://www.ntnu.no).
