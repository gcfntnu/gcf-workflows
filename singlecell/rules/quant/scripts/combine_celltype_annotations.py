#/usr/bin/env python

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse
import logging
import csv

import pandas as pd
import scanpy as sc


def setup_logging(log_file=None):
    handlers = [logging.StreamHandler()]
    if log_file is not None:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers
    )

def load_celltype_annotation(filepath, comment_prefix="#", subset=None):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Detect comments at the top
    start_line = 0
    while start_line < len(lines) and lines[start_line].startswith(comment_prefix):
        start_line += 1

    # Extract a sample of non-commented data for delimiter detection
    sample_data = "".join(lines[start_line : start_line + 5])

    # Detect the delimiter
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(sample_data)
    delimiter = dialect.delimiter

    # Read the file into a pandas DataFrame, skipping detected comments
    df = pd.read_csv(filepath, delimiter=delimiter, skiprows=start_line, index_col=0)
    if df.index.name != "barcode":
        logging.warning(f"Renaming first column `{df.index.name} to `barcode`")
        df.index.name = "barcode"
    if subset is not None:
        df = df[[s for s in subset if s in df.columns]]
    return df
    
def main():
    parser = argparse.ArgumentParser(description="Convert an AnnData h5ad file using a gene mapping file.")
    parser.add_argument("--input", type=str, help="Path to input .h5ad file", required=True)
    parser.add_argument("--annotation", nargs='+', type=str, help="Path to gene mapping tsv file(s)", required=True)
    parser.add_argument("--output", type=str, help="Output .h5ad file", required=True)
    parser.add_argument("--log", default=None, help="Log file")
    args = parser.parse_args()
    
    setup_logging(args.log)
    
    adata = sc.read_h5ad(args.input)

    anno_cols = set()
    for i, anno_fn in enumerate(args.annotation):
        logging.info(f"Procesessing annotations from {anno_fn}")
        anno = load_celltype_annotation(anno_fn)
        anno_cols.update(anno.columns)
        logging.info(f"Adding {anno.shape[1]} annotation columns from {anno_fn}")
        adata.obs = adata.obs.merge(anno, how="left", right_index=True, left_index=True)
    
    logging.info(f"Saving processed data to {args.output}")
    #tiny_data.write_h5ad(args.output, compression='gzip', compression_opts=4)
    merged = adata.obs[list(anno_cols)]
    merged.index.name = "barcode"
    merged.to_csv(args.output, sep="\t")
    logging.info("Processing complete.")

if __name__ == "__main__":
    main()
