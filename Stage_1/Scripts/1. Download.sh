#!/usr/bin/env bash
# Script: download.sh
# Purpose: Download all raw FASTQ files from ENA/SRA using SRR IDs

set -euo pipefail

RAW_DIR="raw"
mkdir -p $RAW_DIR

echo "========================================"
echo " Step: Downloading raw FASTQ files"
echo " Output -> $RAW_DIR"
echo "========================================"

while read SRR; do
  echo ">>> Downloading $SRR ..."
  curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/${SRR: -3}/$SRR/${SRR}_1.fastq.gz \
      -o $RAW_DIR/${SRR}_1.fastq.gz
  curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/${SRR: -3}/$SRR/${SRR}_2.fastq.gz \
      -o $RAW_DIR/${SRR}_2.fastq.gz
done < all_samples.txt

echo ">>> Download finished. Sample of raw files:"
ls -lh $RAW_DIR | head -n 20
