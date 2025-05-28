
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys
import argparse
import logging

import pandas as pd
import scanpy as sc
import anndata as ad

def setup_logging(log_file="script.log"):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def load_gene_map(gene_map_path):
    logging.info(f"Loading gene map from {gene_map_path}")
    return pd.read_csv(gene_map_path, sep="\t", index_col=0)

def load_anndata(h5ad_path):
    logging.info(f"Loading AnnData from {h5ad_path}")
    return sc.read_h5ad(h5ad_path)

def process_features(adata, gene_map, src_organism, dst_organism, non121_strategy="first"):
    logging.info(f"Processing data from {src_organism} to {dst_organism}")
    if src_organism == dst_organism:
        return adata
    merged_var = adata.var.merge(gene_map, left_index=True, right_index=True, how="left")
    merged_var.dropna(axis=0, inplace=True)
    n_gene_id_dups = merged_var.duplicated(subset=f"{dst_organism}_gene_id", keep=False).sum()
    if n_gene_id_dups > 0:
        if non121_strategy == "use_symbol":
            adata = adata[:,merged_var.index]
            adata.var.index = merged_var.index = merged_var[f"{dst_organism}_gene_symbol"]
        else:
            merged_var.drop_duplicates(subset=f"{dst_organism}_gene_id", keep=non121_strategy, inplace=True)
            adata = adata[:,merged_var.index]
            adata.var.index = merged_var.index = merged_var[f"{dst_organism}_gene_id"].copy()
    adata.var = merged_var
    return adata

def compress_data(adata):
    minimal = ad.AnnData(adata.X.astype('f').asformat('csr'),
                         var=pd.DataFrame(index=adata.var.index),
                         obs=pd.DataFrame(index=adata.obs.index)
                         )
    return minimal
    
def main():
    parser = argparse.ArgumentParser(description="Convert an AnnData h5ad file using a gene mapping file.")
    parser.add_argument("--input", type=str, help="Path to input .h5ad file")
    parser.add_argument("--gene-map", type=str, help="Path to gene mapping .tsv file")
    parser.add_argument("--src-organism", type=str, help="Source organism")
    parser.add_argument("--dst-organism", type=str, help="Destination organism")
    parser.add_argument("--output", type=str, help="Output .h5ad file")
    parser.add_argument("--log", type=str, default="mmcell_input.log", help="Log file")
    args = parser.parse_args()
    
    setup_logging()
    
    adata = load_anndata(args.input)
    gene_map = load_gene_map(args.gene_map)
    logging.info(f"Procesessing data ...")
    processed_data = process_features(adata, gene_map, args.src_organism, args.dst_organism)
    tiny_data = compress_data(processed_data)
    logging.info(f"Saving processed data to {args.output}")
    tiny_data.write_h5ad(args.output, compression='gzip', compression_opts=4)
    logging.info("Processing complete.")

if __name__ == "__main__":
    main()
