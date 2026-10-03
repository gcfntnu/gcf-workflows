# single-cell
Single cell RNA seq analysis, 10Xgenomics platform

### Introduction
single-cell: Single cell RNA seq analysis, 10Xgenomics platform

The pipeline is built using [Snakemake](https://bitbucket.org/snakemake/snakemake), a flexible pipeline tool. It runs within docker containers making installation trivial and results highly reproducible.


### Documentation

No documentation and a fair warning that this is still in development and **not** at the monent useful for anybody outside the Genomics Core Facility.

### BFQ deliverables and optional notebooks

For SplitPipe (`quant.method: splitpipe`) and 10x Cell Ranger
(`quant.method: cellranger`), the `bfq_all` and `multiqc_report` targets do not
require the Scanpy preprocessing notebook, its HTML report, or its
`{aggr_id}_preprocessed.h5ad` output. A missing or failed preprocessing notebook
therefore does not block these default BFQ targets.

The canonical `{aggr_id}_filtered.h5ad` produced by `scanpy_aggr_finalize` remains
a required BFQ expression deliverable, along with the existing expression
matrices, QC dependencies, summaries, logs and figures. SplitPipe supplies its
native figures; Cell Ranger's BFQ UMAP is generated independently from the
filtered H5AD by `plotpca.py`.

The notebook rules remain available for optional use. From a prepared project
workdir, using the same Snakemake/container options as the normal workflow run,
request the collector explicitly:

```bash
snakemake --use-singularity --cores 24 bfq_level2_notebooks
```

This collects HTML and executed IPYNB files for every configured aggregate under
`data/tmp/singlecell/bfq/notebooks/` (with the default interim path). To request
only one aggregate's preprocessed H5AD, target its file directly, for example:

```bash
snakemake --use-singularity --cores 24 \
  data/tmp/singlecell/quant/aggregate/splitpipe/scanpy/all_samples_preprocessed.h5ad
```

Use `cellranger` in that path for Cell Ranger, and replace `all_samples` with the
configured aggregate ID. Optional notebook execution still requires compatible
Snakemake and notebook-container versions; removing it from the defaults does
not repair notebook runtime incompatibilities. The separate Parse STARsolo BFQ
defaults are unchanged.

For recovery, retain the existing scientific outputs and Snakemake metadata.
Apply this target-list change to the workflow copy actually used by the project,
then inspect a dry run of `multiqc_report`. A BFQ resume that preserves the
working workflow will not pick up changes made only in `/opt/gcf-workflows`.
Existing optional notebook files are not deleted by changing the default targets.

The focused DAG checks use the real singlecell workflow and temporary copies of
the bundled test inputs. They check both default targets, optional notebook
targets, and retained outputs with a failed notebook log. No scientific jobs or
container downloads run:

```bash
python -m unittest discover -s singlecell/.tests -p 'test_*.py'
```
