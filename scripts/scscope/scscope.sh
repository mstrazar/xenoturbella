#!/usr/bin/env bash

# Run scScope as a script ; red arguments
EXP_ID=$1
BATCH=$2
DATA=$3

# Parameters
DEPTH=128

# Parameters
SCIMP=~/Dev/sc-imputation/code/models/
IN=~/Dev/data/xeno/output/$BATCH/h5ad//$DATA.h5ad

# Make output
OUT_DIR=~/Dev/data/xeno/output/$BATCH/$EXP_ID/h5ad/
mkdir -p $OUT_DIR
echo "Writing outout to $OUT_DIR"

# Start processing
date;hostname;pwd
PY=~/Dev/py3/bin/python
$PY ./scScope.py --input $IN --output-dir $OUT_DIR --batch-var orig.ident --n-latent 50 --depth-encoder $DEPTH --depth-decoder $DEPTH
date