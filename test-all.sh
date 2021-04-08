#!/bin/bash

BASE_DIR=$(pwd)

for PIPELINE in "small-rna" "single-cell" "rna-seq" "microbiome" "wgs"
do
    cd $BASE_DIR/$PIPELINE/.tests/
    snakemake -n -p -j1 bfq_all
    rm -rf .snakemake/
done
