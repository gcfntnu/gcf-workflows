#!/bin/bash

BASE_DIR=$(pwd)

for PIPELINE in "smallrna" "singlecell" "rnaseq" "microbiome" "default" 
do
    cd $BASE_DIR/$PIPELINE/.tests/
    snakemake -n -p -j2 bfq_all
    rm -rf .snakemake/
done
