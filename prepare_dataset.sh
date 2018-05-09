#!/bin/bash
#parameters: read_basedir, reference_fasta, result_dir, SNP_percentage

READ_BASEDIR=$1
REFERENCE_FASTA=$2
RESULT_DIR=$3
SNP_PERCENTAGE=$4

set -e
mkdir $RESULT_DIR

REFERENCE_COPY=$RESULT_DIR/original_reference.fasta
cp $REFERENCE_FASTA $REFERENCE_COPY
bwa index $REFERENCE_COPY 2> /dev/null

for read in $READ_BASEDIR/*.fast5
do
  cp $read $RESULT_DIR/$(basename $read)
done


tombo resquiggle $RESULT_DIR $REFERENCE_COPY --overwrite --bwa-mem-executable bwa --quiet 2> /dev/null
mkdir $RESULT_DIR/reads

for read in $RESULT_DIR/*.fast5
do
  filename=`echo "$read" | sed -r "s/.+\/(.+)\..+/\1/"`
  mkdir $RESULT_DIR/reads/$filename
  mv $read $RESULT_DIR/reads/$filename/read.fast5
done

python3 make_SNPs.py $RESULT_DIR/reads $REFERENCE_COPY $SNP_PERCENTAGE

for read_dir in $RESULT_DIR/reads/*
do
  bwa index $read_dir/reference.fasta 2> /dev/null
  tombo resquiggle $read_dir $read_dir/reference.fasta --overwrite --bwa-mem-executable bwa --quiet 2> /dev/null
done