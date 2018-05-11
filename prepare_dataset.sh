#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "usage: $0 read_basedir reference_fasta result_dir SNP_percentage"
    exit 1
fi

READ_BASEDIR=$1
REFERENCE_FASTA=$2
RESULT_DIR=$3
SNP_PERCENTAGE=$4

set -e
mkdir $RESULT_DIR

echo "[prepare]: initial reference indexing..."
REFERENCE_COPY=$RESULT_DIR/original_reference.fasta
cp $REFERENCE_FASTA $REFERENCE_COPY
bwa index $REFERENCE_COPY #2> /dev/null

echo "[prepare]: copying reads..."

for read in $READ_BASEDIR/*.fast5
do
  cp $read $RESULT_DIR/$(basename $read)
done

echo "[prepare]: initial resquiggle..."
tombo resquiggle $RESULT_DIR $REFERENCE_COPY --overwrite --bwa-mem-executable bwa #--quiet 2> /dev/null

echo "[prepare]: preparing directory for reads..."
mkdir $RESULT_DIR/reads

echo "[prepare]: preparing subdirectories for reads..."
for read in $RESULT_DIR/*.fast5
do
  filename=`echo "$read" | sed -r "s/.+\/(.+)\..+/\1/"`
  mkdir $RESULT_DIR/reads/$filename
  mv $read $RESULT_DIR/reads/$filename/read.fast5
done

echo "[prepare]: making artificial SNPs..."
python3 make_SNPs.py $RESULT_DIR/reads $REFERENCE_COPY $SNP_PERCENTAGE

echo "[prepare]: indexing and resquiggling to changed reference..."
CURRENT_READ=0
NUMBER_OF_READS=`ls $RESULT_DIR/reads -1 | wc -l`
for read_dir in $RESULT_DIR/reads/*
do
  CURRENT_READ=`expr $CURRENT_READ + 1`
  echo "processing read $CURRENT_READ / $NUMBER_OF_READS"
  set +e
  bwa index $read_dir/reference.fasta 2> /dev/null
  set -e
  tombo resquiggle $read_dir $read_dir/reference.fasta --overwrite --bwa-mem-executable bwa --quiet 2> /dev/null
done