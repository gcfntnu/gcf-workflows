# Aggregate AnnData contract

`{aggr_id}_filtered.h5ad` is the canonical quantifier-filtered, cell-called
count dataset. `all_samples_filtered.h5ad` is its all-samples instance. The
aggregate converter retains the quantifier-filtered cells and measured feature
universe, stores raw counts in `X`, and attaches sample, barcode, feature,
doublet, demultiplexing, and other upstream technical metadata when configured.
It does not apply auto-QC, doublet, gene, or analysis selection. For Parse
STARsolo, the established R/T barcode collapse is an upstream identity step:
two technical R/T barcodes can represent one biological cell. The resulting
index uses the canonical Parse barcode and retains the summed raw counts.

Auto-QC produces a separate cell table and pass mask under `auto_qc/`. These
results are not required to build the filtered H5AD. When `qc.qc_sample`
contains `cell_class`, QC explicitly depends on MapMyCells output. General
annotation remains an analysis branch and is not needed to build the filtered
H5AD. Pseudobulk consumes filtered counts and its annotation sidecar directly.

Preprocessing will consume the filtered H5AD and may use the separate QC and
annotation products to create `{aggr_id}_preprocessed.h5ad`. Both H5ADs are
intended as deliverables and independent starting points for downstream work.
The preprocessing rule and analysis choices are not implemented here.

## Aggregate path audit

- Required: quantifier filtered matrices, reference feature metadata, and the
  quantifier barcode sidecar. The Cell Ranger aggregate CSV is required when
  Cell Ranger itself aggregated libraries. Parse STARsolo requires its well/barcode
  sidecar for R/T identity collapse.
- Optional upstream inputs: doublet ranking/classification, multiplex calls,
  CellBender expression presence, CellBender matrix selection, and velocity
  layers. Barcode sidecars merge by canonical barcode; sample identity comes
  from the readers and quantifier sidecars. 10x STARsolo library order is
  validated against the matrix paths before assigning numeric suffixes.
- `convert_scanpy.py` already reads quantifier matrices, aligns and concatenates
  features, merges barcode/feature metadata, adds technical flags, and normalizes
  storage. The retired finalizer uniquely added strict QC alignment and hard
  `autoqc_pass` selection, post-QC CellTypist merging, and Parse R/T collapse.
  QC selection and annotation do not belong in filtered assembly. Parse R/T
  collapse remains in `postprocess_starsolo_rt.py`; the converter creates its
  technical R/T input and the collapse writes the canonical filtered output.
- QC metrics and masks remain separate products. This keeps the filtered
  H5AD independent of the QC DAG, including the optional MapMyCells dependency
  used only when QC stratifies by `cell_class`.
