#!/bin/bash

# 2020-03-8 Wrapper script to run docker container with image containing
# Swift Accel-Amplicon analysis workflow adapted for any reference genome.
# S. Sandhu, S. Chaluvadi, and J. Irish
# Swift Biosciences, Inc.

optflags="$1"
masterfile="$2"
ref_fasta="/refgenomes/covid19.fasta" # path in docker

time docker run -P --rm -e LOCAL_USER_ID=$(id -u $USER) \
    -v ${PWD}:/data \
    -t gcfntnu/sc2:parallel-3 \
    /usr/local/pipelines/swift_amplicon_analysis_sarscov2.sh \
    "$optflags" "$masterfile" "$ref_fasta"

