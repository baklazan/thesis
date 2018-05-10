#!/bin/bash
#parameters: configuration, dataset, output_directory

CONFIG=$1
DATASET=$2
OUTPUT_DIR=$3

set -e
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/results

for dir in $DATASET/reads/*
do
  python3 detect_SNP.py $dir/reference.fasta $dir -i -o $OUTPUT_DIR/results/$(basename $dir).txt -c $CONFIG
done
