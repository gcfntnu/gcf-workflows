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

## Implementation checkpoint

This change implements configuration and selection only. TCR configuration validation occurs before WT-specific rule includes. A valid TCR configuration then stops with `TCR execution rules are not implemented yet`; it cannot run WT quantification accidentally. Dedicated validate/process/combine/link rules and TCR final targets are the next implementation step. Do not deploy this branch for production TCR execution yet.

Run the configuration tests with:

```bash
python -m unittest discover -s singlecell/.tests -p test_immune_config.py -v
```

The existing CI dry-run command is `snakemake -n -p --cores 1 multiqc_report` from `singlecell/.tests`. Use its normal fixture in addition to the explicit TCR stopping check. Configmaker's `LibprepConfig.select` can check the new entry independently without generating a project or inspecting real sequencing data.
