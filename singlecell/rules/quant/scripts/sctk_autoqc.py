import warnings
import gc
import os
import sys
import io
import logging
import argparse
from contextlib import redirect_stdout

import numpy as np
import scanpy as sc
import pandas as pd
from scipy.special import logit, expit
import sctk

warnings.filterwarnings("ignore")

def setup_logger():
    #fixme: this does not vibe very well with capture_prints_to_logger
    logger = logging.getLogger("sctk_autoqc_logger")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    
    print_handler = logging.StreamHandler(sys.stdout)
    print_handler.setLevel(logging.DEBUG)
    print_handler.addFilter(lambda record: record.levelno <= logging.INFO)
    print_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    
    return logger

logger = setup_logger()

def capture_prints_to_logger(func, *args, **kwargs):
    buffer = io.StringIO()
    result = None
    with redirect_stdout(buffer):  # Redirect prints to buffer
        result = func(*args, **kwargs)
    output = buffer.getvalue().strip()
    if output:  # Log the captured output
        for line in output.split("\n"):
            logger.debug(line)
    return result


def inverse_asin(x):
    """
    Compute the inverse arcsin square root transformation.
    
    Parameters
    ----------
    x : array-like
        Transformed values.
    
    Returns
    -------
    array-like
        Original values before transformation.
    """
    return np.square(np.sin(x))

def gcf_qc_sample(adata, qc_sample='sample_id_x_sublib', inplace=True):
    """
    Generate a QC sample ID from a column or by combining two columns in `adata.obs`.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    qc_sample : str, optional
        Column names to combine (separated by "_x_"). Default is "sample_id_x_sublib".
    inplace : bool, optional
        If True, modifies `adata` in place. If False, returns the QC sample ID.
    
    Returns
    -------
    None or pd.Series
        If `inplace` is False, returns the generated QC sample ID.
    """
    
    if "_x_" in qc_sample:
        col1, col2 = qc_sample.split("_x_")
        if col1 in adata.obs.columns and col2 in adata.obs.columns:
            qc_sample_id = adata.obs[col1].astype(str) + "__" + adata.obs[col2].astype(str)
            qc_sample_id = qc_sample_id.astype("category")
        else:
            raise KeyError(f"Both '{col1}' and '{col2}' must be present in adata.obs.")
    elif qc_sample in adata.obs:
        adata.obs["qc_sample_id"] = adata.obs['sample_id'].copy().astype('category')
    else:
        raise KeyError(f"'{qc_sample}' must be present in adata.obs.")
    if inplace:
        adata.obs["qc_sample_id"] = qc_sample_id
    else:
        return qc_sample_id
        
def gcf_calculate_qc_metrics(adata, layers=["counts"], fraction_transforms=["logit"]):
    """
    Compute quality control metrics for an AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    layers : list of str, optional
        Layers to calculate QC metrics for. Defaults to ["counts"].
    fraction_transforms : list of str, optional
        List of fraction transforms to apply ("logit" or "asin").
    """
    pc = adata.var["gene_biotype"].eq("protein_coding") if "gene_biotype" in adata.var.columns else None
    if pc is not None:
        adata.var["pc"] = pc
    qc_vars = [v for v in ["mt", "ribo", "hb", "pc"] if v in adata.var.columns]
    
    for layer in layers:
        layer_key = "" if layer == "counts" else layer
        layer = None if (layer == "counts" and "counts" not in adata.layers) else layer
        obs, _ = sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=False, log1p=True, percent_top=[50], layer=layer)
        
        for v in qc_vars:
            obs[f"{v}_fraction"] = obs[f"pct_counts_{v}"] / 100.0
            
            if "logit" in fraction_transforms:
                obs[f"logit_{v}_fraction"] = logit(np.clip(obs[f"{v}_fraction"], 1e-3, 1 - 1e-3))
            if "asin" in fraction_transforms:
                obs[f"asin_{v}_fraction"] = np.arcsin(np.sqrt(obs[f"{v}_fraction"]))
        
        for col in obs.columns:
            key = f"{layer_key}_{col}" if layer_key else col
            adata.obs[key] = obs[col]
    
    for v in ["background_fraction", "nuclear_fraction"]:
        if v in adata.obs:
            if "logit" in fraction_transforms:
                adata.obs[f"logit_{v}"] = logit(np.clip(adata.obs[v], 1e-3, 1 - 1e-3))
            if "asin" in fraction_transforms:
                adata.obs[f"asin_{v}"] = np.arcsin(np.sqrt(adata.obs[v]))

def extract_normalized_gene_expression(adata, gene="MALAT1"):
    """
    Extract normalized gene expression values for a given gene.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    gene : str, optional
        Gene name or ID to extract. Default is "MALAT1".
    """
    dat = adata.copy()
    sc.pp.normalize_total(dat, target_sum=10000)
    sc.pp.log1p(dat)
    
    if gene in adata.var.index:
        gene_id = gene
    elif "gene_symbols" in dat.var.columns and gene in dat.var["gene_symbols"].values:
        gene_id = dat.var.index[dat.var["gene_symbols"] == gene].values[0]
    else:
        raise ValueError("Gene not found in AnnData object.")
    
    adata.obs[gene] = dat[:, gene_id].X.toarray().squeeze()
    del dat
    gc.collect()

def gcf_prefilter(adata, drop_doublets=True, only_protein_coding=True, update_nuclear_fraction=True, min_genes=100, min_cells=3):
    """
    Filter out low-quality cells and genes based on expression thresholds.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    min_genes : int, optional
        Minimum number of genes expressed per cell. Default is 200.
    min_cells : int, optional
        Minimum number of cells expressing a gene. Default is 3.
    """
    dat = adata.copy()

    if only_protein_coding and 'gene_biotype' in dat.var.columns and 'protein_coding' in dat.var["gene_biotype"]:
        keep_vars = dat.var["gene_biotype"]=="protein_coding"
        logger.debug(f"Subsetting features to protein coding genes ({sum(keep_vars)}/{dat.shape[1]})")
        dat = dat[:, keep_vars]
        if update_nuclear_fraction and 'spliced' in dat.layers and 'unspliced' in dat.layers:
            logger.debug(" Updating nuclear fracion (exon fration of spliced rna) based on only protein coding genes")
            u_sum = np.asarray(dat.layers['unspliced'].sum(1)).squeeze()
            s_sum = np.asarray(dat.layers['spliced'].sum(1)).squeeze()
            dat.obs['nuclear_fraction'] = u_sum / (u_sum + s_sum) # re-calculate nuclear fraction based on protein coding genes
    if drop_doublets and "droplet_type" in dat.obs.columns:
        keep_obs = dat.obs["droplet_type"]=="singlet"
        logger.debug(f"Subsetting barcodes to singlets only ({sum(keep_obs)}/{dat.shape[1]})")
        dat = dat[dat.obs["droplet_type"]=="singlet"]
    n_obs, n_vars = dat.shape
    sc.pp.filter_cells(dat, min_genes=min_genes)
    if n_obs != dat.shape[0]:
        logger.debug(f"Subsetting barcodes based on min_genes requirement ({dat.shape[0]}/{n_obs})")
    sc.pp.filter_genes(dat, min_cells=min_cells)
    if n_vars != dat.shape[1]:
        logger.debug(f"Subsetting genes based on min_cells requirement ({dat.shape[1]}/{n_vars})")
    prefilter_mask = adata.obs.index.isin(dat.obs.index) # 1: keep 0: remove
    logger.debug(f"Total cells after prefilter: ({dat.shape[0]}/{adata.shape[0]})")
    logger.debug(f"Total genes after prefilter: ({dat.shape[1]}/{adata.shape[1]})")
    adata.obs["prefilter_mask"] = prefilter_mask
    logger.debug(adata.obs.groupby('qc_sample_id')['prefilter_mask'].sum().to_string())
    logger.debug(adata.obs.groupby('qc_sample_id')['droplet_type'].value_counts().to_string())
    
def gcf_cellwise_qc(adata, metrics, cell_qc_key="cell_passed_qc", uns_qc_key="scautoqc_ranges", **kwargs):
    """
    Compute cell-wise QC metrics and add them to `adata.obs`.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    qc_vars : list of str
        List of QC variables to compute.
    """
    from sctk._pipeline import fit_gaussian
    
    metric_params = metrics.loc[:, ["min", "max", "scale", "side", "min_pass_rate"]].replace({np.nan: None}).T.to_dict(orient="list")
    n_obs = adata.n_obs
    range_df = pd.DataFrame(columns=["low", "high"])
    pass_filter = {}
    
    for m, v in metric_params.items():
        min_x, max_x, scale, side, min_pass_rate = v
        if m not in adata.obs.columns:
            continue
        x = adata.obs[m].values.astype(np.float32)
        if scale == "log":
            x, min_x, max_x = np.log1p(x), np.log1p(min_x) if min_x else None, np.log1p(max_x) if max_x else None
        elif scale == "logit":
            x, min_x, max_x = logit(x), logit(min_x) if min_x else None, logit(max_x) if max_x else None
        elif scale == "asin":
            x, min_x, max_x = np.arcsin(np.sqrt(x)), np.arcsin(np.sqrt(min_x)) if min_x else None, np.arcsin(np.sqrt(max_x)) if max_x else None
        
        try:
            x_low, x_high, _ = capture_prints_to_logger(fit_gaussian, x, xmin=min_x, xmax=max_x, **kwargs)
        except ValueError:
            logging.debug("--::-- ")
            x_low, x_high = min_x if min_x else x.min(), max_x if max_x else x.max()
        
        pass_filter[m] = (x_low <= x) & (x <= x_high)
        range_df.loc[m] = [expit(x_low) if scale == "logit" else inverse_asin(x_low) if scale == "asin" else np.expm1(x_low) if scale == "log" else x_low,
                            expit(x_high) if scale == "logit" else inverse_asin(x_high) if scale == "asin" else np.expm1(x_high) if scale == "log" else x_high]
        logger.info(f"cutoffs qc var: {m}    low={x_low} high={x_high}")
    
    all_passed = np.all(list(pass_filter.values()), axis=0)
    adata.obs[cell_qc_key] = all_passed
    adata.uns[uns_qc_key] = range_df
    logger.info(range_df)

def gcf_default_metrics(qc_vars):
    """
    Define default QC metric thresholds.
    
    Returns
    -------
    metrics: pd.DataFrame
        Dataframe of default QC thresholds.   
    """
    metrics = pd.DataFrame([(200, None, "log", "min_only", 0.1),
                            (100, None, "log", "min_only", 0.1),
                            (0, 0.2, "logit", "min_only", 0.1),
                            (0.2, 0.9, "logit", "both", 0.1)],
                           index = ["total_counts", "n_genes_by_counts", "mt_fraction", "nuclear_fraction"],
                           columns = ["min", "max", "scale", "side", "min_pass_rate"])
    if len(set(qc_vars).intersection(metrics.index.values)) < len(qc_vars):
        uniq_qc_vars = list(set(qc_vars).intersection(metrics.index.values))
        logger.critical(f" {str(uniq_qc_vars)} are named qc variables not present in default metrics")
        raise ValueError("There are qc_vars element not valid for default metrics")
    
    return metrics.loc[qc_vars,:]

def gcf_run_autoqc(adata, qc_vars, qc_sample='sample_id_x_sublib'):
    """
    Run automatic QC filtering on an AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    qc_vars : list of str
        List of QC variables to use for filtering.
    qc_sample: str
        QC Sample ID identifier
    """
    gcf_qc_sample(adata, qc_sample=qc_sample)
    metrics = gcf_default_metrics(qc_vars)
    metrics_list = metrics["scale"] + "_" + metrics.index
    metrics_list = [str(i).replace("log_", "log1p_") for i in metrics_list]
    gcf_prefilter(adata)
    gcf_calculate_qc_metrics(adata)
    consensus_passed_qc = pd.Series([False]*adata.shape[0], index=adata.obs.index, name="sctk_autoqc_mask")
    consensus_passed_qc.index.name = "barcode"
    for qc_sample in adata.obs['qc_sample_id'].unique():
        cell_subset =  (adata.obs['qc_sample_id']==qc_sample) & (adata.obs["prefilter_mask"])
        if cell_subset.sum() < 100:
            logger.debug(f"Extremely small number of valid barcodes in qc_sample_id: {qc_sample}: {cell_subset.sum()}")
            logger.debug(f"Number of QC Sample barcodes: {sum(adata.obs['qc_sample_id']==qc_sample)}")
            logger.debug(f"Number of barcodes after prefilter: {sum(adata.obs['prefilter_mask'])}")
                  
        dat = adata[cell_subset,:]
        gcf_cellwise_qc(dat, metrics=metrics)
        #capture_prints_to_logger( sctk.multi_resolution_cluster_qc, dat, metrics=metrics_list )
        sctk.multi_resolution_cluster_qc(dat, metrics=metrics_list)
        
        passed = dat.obs['consensus_passed_qc']
        ok_cells = passed.index[passed].values
        logger.info(f" PASSED {qc_sample}: {sum(passed)}/{dat.shape[0]}")
        consensus_passed_qc[ok_cells] = True
    logger.info(f" PASSED TOTAL: {sum(consensus_passed_qc)}/{len(consensus_passed_qc)}")
    return consensus_passed_qc

def check_args(args):
    if args.qc_sample is None or args.qc_sample == 'auto':
        if args.quantifier.startswith('parsebio') or args.quantifier.startswith('splitpipe'):
            args.qc_sample = 'sample_id_x_sublib'
        elif args.quantifier.startswith('cellranger'):
            args.qc_sample = 'sample_id_x_flowcell_id'
        else:
           args.qc_sample = 'sample_id' 
    
    args.qc_vars = args.qc_vars.split(',')
    return args


def qc_plotter(adata):
    pass

def main():
    valid_quantifiers = ['parsebio_starsolo', 'parsebio_starsolo_cellbender', 'splitpipe', 'splitpipe_cellbender', 'cellranger', 'cellranger_cellbender']
    parser = argparse.ArgumentParser(description="Perform quality control on an AnnData object.")
    parser.add_argument("-i", "--input", type=str, help="Path to input AnnData (.h5ad) file")
    parser.add_argument("-o", "--output", type=str, help="Path to save the output AnnData (.h5ad) file")
    parser.add_argument("--quantifier", type=str, default="parsebio_starsolo", help="Gene count quantifier method", choices=valid_quantifiers)
    parser.add_argument("--qc-sample", type=str, default="sample_id_x_sublib", help="Division of barcodes where QC is run independently")
    parser.add_argument("--qc-vars", type=str, default="total_counts,n_genes_by_counts,mt_fraction,nuclear_fraction", help="QC variables")
    parser.add_argument("--log-filename", type=str, default="sctk_autoqc.log", help="Log file")
    args = parser.parse_args()
    args = check_args(args)

    if args.log_filename:
        log_handler = logging.FileHandler(os.path.join(os.path.dirname(args.output), args.log_filename), mode="w")
        log_handler.setLevel(logging.DEBUG)
        log_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(log_handler)

    logger.info(" Loading input AnnData file...")
    adata = sc.read_h5ad(args.input)
    adata.X = adata.X.astype('f')

    logger.info(" Performing QC computations...")
    passed = gcf_run_autoqc(adata, qc_vars=args.qc_vars, qc_sample=args.qc_sample)

    logger.info(" Saving passed barcode mask ...")
    passed.astype(int).to_csv(args.output, float_format="%d", sep="\t")
    assert(passed.shape[0] == adata.shape[0])

if __name__ == "__main__":
    main()
    logger.info(" QC process completed successfully.")
