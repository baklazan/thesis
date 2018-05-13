#!/bin/bash
#parameters: configuration, dataset, output_directory

if [ "$#" -ne 3 ]; then
    echo "usage: $0 configuration dataset output_directory"
    exit 1
fi

CONFIG=$1
DATASET=$2
OUTPUT_DIR=$3

set -e
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/results

echo -n $DATASET > $OUTPUT_DIR/dataset.txt

set +e
for dir in $DATASET/reads/*
do
  python3 detect_SNP.py $dir/reference.fasta $dir -o $OUTPUT_DIR/results/$(basename $dir).txt -c $CONFIG
done
