#!/bin/bash

BASE_DIR=$(pwd)

for PIPELINE in "smallrna" "singlecell" "metagenome" "rnaseq" "microbiome" "default" "wgs" 
do
    cd $BASE_DIR/$PIPELINE/.tests/
    snakemake -n -p -j2 bfq_all
    rm -rf .snakemake/
done
