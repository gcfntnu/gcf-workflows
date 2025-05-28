#/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import os
import argparse
import logging
import csv
import pickle

import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
try:
    import rapids_singlecell as rsc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=False,  # Allows oversubscription
        pool_allocator=False,  # default is False
        devices=0,  # GPU device IDs to register. By default registers only GPU 0.
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
except ImportError:
    pass

def setup_logging(log_file="script.log"):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def guess_org_from_model(model):
    if 'human' in model.lower():
        return 'homo_sapiens'
    elif 'mouse' in model.lower():
        return 'mus_musclus'
    elif 'bts_atlas' in model.lower():
        return 'homo_sapiens'
    elif 'immune' in model.lower():
        return 'homo_sapiens'
    else:
        with open(model, "rb") as fh:
            try:
                pkl_obj = pickle.load(fh)
                desc = pkl_obj['description']
            except:
                pass
            guess_org_from_model(desc)
    return None

def run_normalize_and_annotate(adata, args):
    # normalize data
    adata_copy = adata.copy()
    logging.info(f"Log transforming and normalizing data ...\n")
    sc.pp.filter_genes(adata_copy, min_cells=3)
    sc.pp.normalize_total(adata_copy, target_sum = 1e4)
    sc.pp.log1p(adata_copy)
    
    # run annotation
    result = celltypist.annotate(adata_copy, model=args.model, majority_voting=True, use_GPU=False)
    #lr_classifier = models.Model.load(args.model)
    #clf = celltypist.classifier.Classifier(adata_copy, lr_classifier)
    #result = clf.celltype(mode = args.mode.replace('_', ' '), p_thres = 0.5)
    #result.predicted_labels['conf_score'] = result.probability_matrix.max(axis=1).values
    return result

def batch_majority_vote(adata, pred, args):
    #adata_copy = adata.copy()
    if args.use_GPU:
        logging.info(f"Using gpu on rapids: v {rsc.__version__}")
        logging.info(f"Log transforming and normalizing data on GPU ...\n")
        adata.X = adata.X.astype('f')
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"])
        rsc.pp.filter_genes(adata, qc_var='n_cells_by_count', min_count=3)
        rsc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', batch_key=args.batch)
        adata = adata[:, adata.var.highly_variable]
        rsc.pp.normalize_total(adata, target_sum = 10000)
        rsc.pp.log1p(adata)
        #rsc.pp.scale(adata, max_value=10)
        rsc.pp.pca(adata, n_comps=50)
        rsc.pp.harmony_integrate(adata, key=args.batch)
        rsc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
        rsc.get.anndata_to_CPU(adata, convert_all=True)
    else:
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum = 1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes = 2000, flavor='seurat_v3', batch_key=args.batch)
        #sc.pp.scale(adata, max_value=10)
        sc.pp.pca(adata, n_comps=50)
        sc.external.pp.harmony_integrate(adata, key=args.batch)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")

        
    lr_classifier = models.Model.load(args.model)
    clf = celltypist.classifier.Classifier(adata, lr_classifier)
    clusters = clf.over_cluster(use_GPU=args.use_GPU)
    return clf.majority_vote(pred, clusters)
        
def main():
    parser = argparse.ArgumentParser(description="Convert an AnnData h5ad file using a gene mapping file.")
    parser.add_argument("--input", type=str,
                        help="Path to input .h5ad file", required=True)
    parser.add_argument("--output", type=str,
                        help="Output tsv file", required=True)
    parser.add_argument("--model", type=str,
                        help="Celltypist model (.pkl)", required=True)
    parser.add_argument("--qc-mask", type=str,
                        help="Barcode QC mask")
    parser.add_argument("--gene-map", type=str, 
                        help="Path to gene mapping .tsv file")
    parser.add_argument("--src-organism", type=str, default="homo_sapiens",
                        help="Source organism")
    parser.add_argument("--dst-organism", type=str, default="homo_sapiens",
                        help="Destination organism")
    parser.add_argument("--batch", default=None,
                        help="Name of column in adata.obs containg the batch variable. Will run annotation separately within each batch.\
                        The majority voting refinement (overclustering) is performed on the whole dataset batch corrected with harmonypy")
    parser.add_argument("--mode", default="best_match",  choices=['best_match','prob_match'],
                        help="Celltypist's mode arg")
    parser.add_argument("--use-GPU", action="store_true", default=False,
                        help="Whether to use GPU for over clustering on the basis of `rapids-singlecell`")
    parser.add_argument("--plot", action="store_true", default=False,
                        help="Whether to save figures from celltypist")
    parser.add_argument("--log", type=str, default="auto_annotate_scanpy.log",
                        help="Log file")
    parser.add_argument("--no-qc-filter", action="store_true", default=False,
                        help="Activate skip qc-filter (override --qc-mask)")
    args = parser.parse_args()
    if args.use_GPU and 'rapids_singlecell' not in sys.modules:
        logging.warn("Warning: rapids_singlecell is not installed but required for GPU running, will switch back to CPU")
        args.use_GPU = False
        rsc = sc
    setup_logging(args.log)
    
    adata = sc.read_h5ad(args.input)
    original_obs_index = adata.obs.index.copy()
    
    if args.gene_map:
        gene_map = pd.read_table(args.gene_map, index_col=0)
        
    if args.src_organism == args.dst_organism:
        adata.var.index = adata.var["gene_symbol"].astype(str)
    else:
        adata.var = adata.var.merge(gene_map, how="left", right_index=True, left_index=True)
        keep_genes = adata.var[f"{args.dst_organism}_gene_symbol"].notna()
        logging.info(f"{sum(keep_genes)}/{adata.shape[1]} genes remaining after converting to ortholog gene symbol")
        adata = adata[:,keep_genes]
        adata.var.index = adata.var[f"{args.dst_organism}_gene_symbol"].astype(str)
        adata.var_names_make_unique() # we are just ignoring the annotation mismatch from this
    if args.qc_mask and not args.no_qc_filter:
        cell_subset = pd.read_table(args.qc_mask, index_col=0).astype('bool')
        barcodes = cell_subset.loc[cell_subset.values,:].index
        n_keep, n_tot = len(barcodes), int(adata.shape[0])
        logging.info(f" {n_keep} of {n_tot} cells remaining after auto_qc")
        adata = adata[barcodes,:]

    if args.batch is not None:
        try:
            batch_var = adata.obs[args.batch]
        except:
            logging.error(f"Batch variable '{str(args.batch)}' not found in adata.obs.columns")
            raise KeyError
        for i, b in enumerate(batch_var.unique()):
            subset = adata[batch_var==b,:]
            logging.info(f"Running batch '{args.batch}=={b}'. Number of cells: {subset.shape[0]}")
            sub_res = run_normalize_and_annotate(subset, args)
            if i == 0:
                pred = sub_res.predicted_labels
                prob = sub_res.probability_matrix
            else:
                pred = pd.concat([pred, sub_res.predicted_labels], axis='index')
                prob = pd.concat([prob, sub_res.probability_matrix], axis='index')
        sub_res.predicted_labels = pred #lets just monkey patch the last sub result with the concatenated results
        result = batch_majority_vote(adata, sub_res, args)
    else:
        result = run_normalize_and_annotate(adata, args)

    pred = result.predicted_labels
    prob = result.probability_matrix
    pred['majority_voting_conf_score'] = [row[pred.majority_voting[index]] if pred.majority_voting[index] in row.index else row.max() for index, row in prob.iterrows()]
    pred['celltypist_cell_type'] = pred['majority_voting'].copy()
    # expand predictions to input adata and fill missing values sensibly
    all_cells = pd.DataFrame([], index=original_obs_index)
    pred = pred.merge(all_cells, how="right", left_index=True, right_index=True)
    for col in ['celltypist_cell_type', 'predicted_labels']:
        if col in pred.columns:
            pred[col] = pred[col].copy().astype("string").fillna('qc_fail').astype('category')
    if 'over_clustering' in pred.columns:
        pred['over_clustering'] = pred["over_clustering"].copy().astype("string").fillna("-1").astype('category')
    if 'conf_score' in pred.columns:
        pred['conf_score'].fillna(0, inplace=True)
    if 'majority_voting_conf_score' in pred.columns:
        pred['majority_voting_conf_score'].fillna(0, inplace=True)
    logging.info(f"Saving processed data to {args.output}")
    pred.to_csv(args.output, sep="\t")
    
    if args.plot:
        result.to_plots(os.path.dirname(args.output))
    
    logging.info("\nProcessing complete.")

if __name__ == "__main__":
    main()
