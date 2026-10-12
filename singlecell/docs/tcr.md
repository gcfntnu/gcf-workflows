# Parent-linked Parse TCR configuration

The `Parse Biosciences Evercode TCR v4 PE` libprep entry describes human TCR-enriched paired-end libraries associated with an existing WT-size v4 experiment. Select it with configmaker's `--libkit 'Parse Biosciences Evercode TCR v4'`. Actual read lengths come from run metadata; the preset neither trims nor imposes the WT run's read lengths. Mini, Mega, other chemistries and standalone TCR are not supported by this initial contract.

Split-Pipe uses `quant.kit: WT` together with `--tcr_analysis`, not a kit value of TCR. `quant.assay: tcr` selects the intended workflow branch; the default `wt` preserves existing expression workflows. The new `docker.splitpipe_immune` entry pins the tested immune image independently of the existing WT image.

## Project configuration

Merge a project overlay after configmaker's generated settings, before loading workflow rules. Snakemake supports ordered config files, for example `--configfile config.yaml tcr-overlay.yaml`; validate the effective merge during workdir setup. The overlay must supply:

```yaml
immune:
  species: human
  reference: default  # imgt is an explicit alternative
  wt_parent:
    root: /absolute/path/to/completed/wt/quant
    sublibraries:
      sublib1: splitpipe/sublib1
      sublib2: splitpipe/sublib2
    combined: aggregate/splitpipe
  sublibrary_order: [sublib1, sublib2]
```

Every key must match a TCR `samples` entry. Each relative parent path is unique and stays beneath the root. The combined parent path is distinct from individual parents. Explicit order must match the order used by the original WT combine; the structural validator checks completeness and uniqueness, while comparison with WT metadata belongs in the later server-side preflight. All configured paths must be visible inside the actual execution environment, with WT inputs bound read-only.

`quant.aggregate.skip: true` disables generic WT expression aggregation; it does not request skipping the future TCR combine rule. Velocity and CellBender are disabled. Retain all WT sample identities and use sublibrary plus barcode when associating immune results with cells.

The preset deliberately has no project-specific WT paths. The validator does not access the filesystem, infer pairing, validate FASTQs, or assert that biological kit metadata is correct; those checks belong to workdir preflight. Missing mappings, incompatible settings and unsupported assay values fail explicitly.

## Workflow targets and execution

The TCR branch includes only immune rules; the WT branch retains its existing includes.
Always pass an explicit target:

- `tcr_preflight_all`: validate the external WT parents, verify every FASTQ MD5 and every paired record, and run Split-Pipe setup-only checks. This target cannot reach analysis rules.
- `data/tmp/singlecell/immune/sublibraries/sublib1/complete.json`: pilot, using the appropriate interim root for the deployment.
- `tcr_all` or `bfq_all`: process all sublibraries, combine, then create a WT-preserving cell table and overlap report.

The parent validator checks version 1.8.2, kit WT, chemistry v4, consistent sample/well definitions and exact agreement between the ordered individual cell calls and the combined WT cell list. Parent directories are read-only external resources. All metadata consumed by Split-Pipe is declared as input. Geometry/kit assumptions still require laboratory confirmation.

Each run has a separate output directory, command/image/revision/exit provenance, log and benchmark. Existing nonempty output directories are rejected; preserve a failed attempt under a different name before retrying. Snakemake-created empty directories are allowed. Helpers require the installed Split-Pipe 1.8.2 output contract, including compressed AIRR/contig files, barcode and clonotype tables, and an HTML report. Empty-repertoire edge cases must be checked during the real pilot; no real TCR analysis has yet validated those output assumptions.

Combine receives the **combined WT parent** and explicit ordered TCR sublibrary paths. Its setup-only check runs immediately before actual combine, since it requires completed immune outputs. The linked table retains every WT cell, summarizes productive TRA/TRB calls, and leaves the original chain-level AIRR outputs intact. Filtered TCR barcodes absent from WT cause a failure with a saved overlap report. No downstream WT matrix is modified.

Configure scheduler-wide `mem_mb` and `io_jobs` budgets as well as cores. Per-rule resource declarations alone do not enforce OS memory limits; the project launcher also caps its Docker container. Start production with one TCR job until the pilot establishes actual resource use.

Tests:

```bash
python -m unittest discover -s singlecell/.tests -p 'test_*.py' -v
```

The legacy CI command remains `snakemake -n -p --cores 1 multiqc_report` from `singlecell/.tests`. The TCR project supplies a small synthetic parent/FASTQ fixture to check the complete DAG without real sequencing data. FASTQ parser tests require pysam, supplied by the pinned immune image.
