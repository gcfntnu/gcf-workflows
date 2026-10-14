# Parent-linked Parse TCR: reproducible execution

This branch processes human Evercode v4 TCR-enriched paired-end libraries with an existing WT-size whole-transcriptome analysis.
It runs Split-Pipe per sublibrary, combines in the original WT order, and links receptor calls to the existing WT cells.
It does not rerun WT expression analysis or modify WT results.
Mini, Mega, other chemistries and standalone TCR are outside this implementation's supported configuration.

## Required inputs

- A demultiplexed TCR run with paired FASTQs and the usual configmaker run metadata, sample submission and checksum manifest.
  Configmaker reads that metadata; the TCR workflow does not recalculate FASTQ checksums or scan read pairing.
- One completed **Split-Pipe 1.8.2, kit WT, chemistry v4** WT output directory for each TCR sublibrary.
  The TCR-to-WT mapping must reflect the physical libraries; matching names alone cannot establish this.
- The completed combined WT output from exactly those sublibraries, in the recorded original order.
- Preserve the complete WT parent directories as read-only resources, including reports.
  The validator reads `process/run_proc_def.json` and each declared sample's `DGE_filtered/cell_metadata.csv`, including `all-sample`.
  Split-Pipe also uses parent resources when creating its reports; do not construct parent stubs containing only the validator inputs.
  WT FASTQs and BAMs are not direct TCR workflow inputs.

Individual parents must have consistent biological sample/well definitions.
The combined WT cell list must be the exact ordered union of individual lists, including `__s1`, `__s2`, etc. suffixes and sample assignments.
Biological sample names/wells are inherited from WT; `config.samples` instead identifies the **technical TCR sublibraries**.
Do not relabel WT samples or guess missing library relationships to make validation pass.

## Software and reference

`docker.config` pins `docker.splitpipe_immune` to:

```text
gcfntnu/psbio-immune@sha256:53ef5bcd84bc7edf9d8173f48f896c819ce298efb9e87383994a30d77c2cfdbe
```

This is the tested `1.8.2-gcf1` image with Split-Pipe 1.8.2 and IgBlast 1.22.0 and installed immune resources.
Use `immune.reference: default` for the tested human TCR database selection.
`imgt` explicitly adds `--use_imgt_db`; that alternative was not used in the completed run.
No separate WT genome argument is needed: `--parent_dir` supplies the WT experiment settings.
The original WT image entry is unchanged.

The tested orchestration environment contains configmaker 0.3.0, Snakemake 9.7.1 and Apptainer:

```text
gcfntnu/bfq@sha256:cdb9d25fec222b75ca32e9cdbf471a6f3965bcd1986813fa7e3e8bcc7c6b689f
```

## 1. Generate a fresh workdir

Use a new project workdir and deploy this repository at a recorded commit under `src/gcf-workflows`.
The following commands assume configmaker is available in the execution environment and all paths are visible there.
Replace the example paths, project ID and revision before executing; do not regenerate an existing analysis workdir.

```bash
work=/absolute/path/to/project/singlecell
run=/absolute/path/to/demultiplexed/flowcell
revision='REPLACE_WITH_FULL_WORKFLOW_COMMIT'
mkdir -p "$work/src"
cd "$work"
git clone git@github.com:gcfntnu/gcf-workflows.git src/gcf-workflows
git -C src/gcf-workflows checkout --detach "$revision"
git -C src/gcf-workflows rev-parse HEAD > WORKFLOW_COMMIT

configmaker.py --project-id GCF-YYYY-NNN --organism homo_sapiens \
  --libkit 'Parse Biosciences Evercode TCR v4' \
  --libprep-config src/gcf-workflows/libprep.config \
  --libprep-entry 'Parse Biosciences Evercode TCR v4 PE' \
  --output config.yaml "$run"
```

The CLI kit name is unsuffixed; the libprep entry ends in `PE`.
Read lengths come from run metadata; the preset does not impose WT read lengths or trim reads.
The completed TCR run used `--expected-read-geometry 242 58` as a metadata assertion; add it only when appropriate for the new run.
Configmaker generates the Snakefile, PEP/configuration and raw-input links.
Verify that its Snakefile includes `src/gcf-workflows/singlecell/singlecell.smk` and that `config.samples` contains the intended sublibraries.
If the submission has a generic kit label, the explicit libprep entry above selects the TCR preset.

## 2. Add the project overlay

Save this as `tcr-overlay.yaml` next to `config.yaml`, expanding the mapping to include every TCR sublibrary:

```yaml
quant:
  method: splitpipe
  assay: tcr
  kit: WT
  chemistry: v4
  use_velo: false
  aggregate:
    skip: true
  cellbender:
    enabled: false
immune:
  species: human
  reference: default
  wt_parent:
    root: /absolute/path/to/completed/wt/quant
    sublibraries:
      sublib1: splitpipe/sublib1
      sublib2: splitpipe/sublib2
    combined: aggregate/splitpipe
  sublibrary_order: [sublib1, sublib2]
```

Every mapping key must match a `config.samples` key exactly.
Paths below `root` must be relative, distinct and remain beneath that root.
`immune.sublibrary_order` must include every sublibrary exactly once and reproduce the WT combine order.
The combined parent is a separate resource, never an individual WT parent reused for all libraries.

Pass **both** config files on every invocation, in this order: `--configfile config.yaml tcr-overlay.yaml`.
The overlay is intentionally separate from generated configuration and is not loaded automatically.
The explicit settings above agree with the libprep preset and make the effective assay unambiguous.
`quant.kit: WT` plus Split-Pipe's `--tcr_analysis` is correct; do not set the kit to `TCR`.
`quant.aggregate.skip: true` skips generic WT aggregation, **not** TCR combine.
Velocity and CellBender remain disabled. The default `quant.assay: wt` preserves existing expression workflows.

## 3. Configure the runtime and inspect the DAG

Keep project work, FASTQs and WT parents visible at the same absolute paths inside the execution environment.
Bind the project read-write and the FASTQ/WT roots read-only, including targets of any symlinks.
For Apptainer, a native-host example is:

```bash
export GCF_SECRET=""
export GCF_INTERIM="$work/data/tmp"
export GCF_EXT="$work/data/ext"
runtime="$(dirname "$work")/runtime"
wt=/absolute/path/to/completed/wt/quant
mkdir -p "$runtime"/{cache,tmp,images} "$work/logs"
export APPTAINER_CACHEDIR="$runtime/cache"
export APPTAINER_TMPDIR="$runtime/tmp"
export TMPDIR="$runtime/tmp"
export APPTAINER_BINDPATH="$(dirname "$work"):$(dirname "$work"),$run:$run:ro,$wt:$wt:ro"

# Bash array: reuse these arguments for every target.
smk=(snakemake --cores 24 --resources mem_mb=192000 io_jobs=4
     --use-apptainer --apptainer-prefix "$runtime/images"
     --configfile config.yaml tcr-overlay.yaml
     --config "workflow_commit=$(cat WORKFLOW_COMMIT)")
"${smk[@]}" --dry-run --printshellcmds -- tcr_all
```

The `--` terminates Snakemake's variable-length configuration arguments before the target.
`GCF_SECRET` may be empty for this TCR-only branch but must exist for shared configuration loading.
Use an absolute interim path; if configuration explicitly sets `interim_dir`, it takes precedence over `GCF_INTERIM`.
Keep image caches outside analysis output directories.
When running the BFQ image as the outer Docker runtime, make the same mounts in Docker **and** Apptainer,
provide a writable HOME, and mount `/etc/passwd` and `/etc/group` read-only for the invoking numeric UID.
The tested Docker launcher used `--privileged` for nested Apptainer and an outer limit of 24 CPUs/220 GiB RAM.
Set equivalent `SINGULARITY_*` variables when the runtime uses those compatibility names.

## 4. Execute explicit targets

Run on the execution server through its scheduler or in a persistent byobu/tmux session.
At GCF, use the existing `geiramh` session for project jobs and retain stdout/stderr and exit status.

```bash
# Existing parent checks and Split-Pipe setup-only validation; no analysis target is reachable.
"${smk[@]}" -- tcr_preflight_all

# Optional first-sublibrary pilot (adjust when INTERIM_DIR differs).
"${smk[@]}" -- "$GCF_INTERIM/singlecell/immune/sublibraries/sublib1/complete.json"

# All sublibraries, ordered combine and WT linkage. Completed dependencies are reused.
set -o pipefail
"${smk[@]}" --printshellcmds --keep-going --keep-incomplete --latency-wait 60 -- tcr_all \
  2>&1 | tee "logs/tcr-production-$(date -u +%Y%m%dT%H%M%SZ).log"
status=${PIPESTATUS[0]}
printf '%s\n' "$status" > logs/tcr-production.exit-status
```

`tcr_all` and `bfq_all` are equivalent in this branch. Production includes setup-only dependencies automatically;
running `tcr_preflight_all` separately is optional.
Setup-only checks retain Split-Pipe's own chemistry sampling and can take several minutes.
There are no workflow FASTQ checksum, full-record integrity or read-pair verification rules.
Input errors are left to the analysis software to report; no additional preflight checks are introduced.
Existing checks validate configuration, WT metadata, successful completion and expected output contracts.

Each Split-Pipe job requests 12 threads, 96,000 MB and one `io_jobs` slot.
The example budget permits two concurrent jobs; use a smaller budget where needed.
Resource declarations constrain scheduling, not OS memory use.
The completed eight-sublibrary run took approximately 4.5 hours with two concurrent jobs;
temporary disk usage reached roughly 450 GiB. Allow sufficient scratch space for that workload and final outputs.

## Outputs, provenance and retries

Under `$GCF_INTERIM/singlecell/immune` (or the configured `INTERIM_DIR`):

| Location | Contents |
| --- | --- |
| `checks/wt_parents.json` | Effective settings, parent identities, metadata hashes and cell counts |
| `preflight/<sublib>/` | Split-Pipe setup-only output |
| `sublibraries/<sublib>/` | Individual TCR analyses and `complete.json` |
| `combined/` | Ordered immune combine, using the **combined WT parent** |
| `combined/all-sample/TCR_filtered/` | `tcr_annotation_airr.tsv.gz`, `barcode_report.tsv`, `clonotype_frequency.tsv`, `tcr_contigs.fa.gz` |
| `combined/all-sample/TCR_unfiltered/` | Corresponding unfiltered receptor outputs |
| `combined/all-sample_analysis_summary.html` | Split-Pipe combined report |
| `linked/wt_tcr_cells.tsv.gz` | Every WT cell plus TCR contig/productive TRA/TRB counts and paired-chain indicator |
| `linked/overlap.json` | WT/TCR overlap counts and any unmatched filtered TCR cells |
| `logs/`, `benchmarks/`, `provenance/` | Per-step logs, resource measurements, commands, image, workflow commit and exit status |

The linked table retains WT cells without TCR calls and preserves original biological sample labels.
Filtered TCR cells outside the combined WT parent cause failure with an overlap report.
The original AIRR chain records and WT expression matrices remain unchanged.

The runner refuses existing nonempty output directories to prevent accidental mixing of attempts.
For a failed job, inspect its log/provenance, preserve its failed output directory under a different name,
and then rerun the same target with `--rerun-incomplete` before `--`.
Keep the configuration and recorded workflow revision fixed during a run;
changing them may cause Snakemake to rerun existing outputs, which the overwrite guard will reject.
Do not delete successful outputs or force a completed production run merely to test documentation changes.

## Validation history and developer checks

The precursor revision `393f3d631219c91c4ba793affe376553f871032d` completed all eight sublibraries,
combine and WT linkage on 2026-10-12 using the pinned image and default database.
It retained 116,102 WT cells and linked 59,692 filtered TCR cells with zero unmatched cells.
The final integration removes the unused optional FASTQ integrity scanner and documents the procedure;
analysis commands, resource requests and output contracts are preserved.
This establishes the tested parent-linked path, not coverage of every empty-repertoire or alternate-reference case.

Run configuration, parent-linkage and runner regression tests from the repository root:

```bash
python -m unittest discover -s singlecell/.tests -p 'test_*.py' -v
```

The legacy WT DAG check remains `snakemake -n -p --cores 1 multiqc_report` from `singlecell/.tests`
with the normal GCF environment configured. Inspect a TCR DAG with both config files before launching a new project.
