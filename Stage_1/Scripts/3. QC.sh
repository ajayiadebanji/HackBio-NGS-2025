#!/bin/bash
# Script: 03_qc.sh
# Purpose: Run FastQC + MultiQC on selected reads
# Why:    Detect sequencing quality issues before trimming.

SELECTED_DIR="selected"
QC_DIR="qc/raw_fastqc"
MULTIQC_DIR="qc/multiqc_raw"

mkdir -p $QC_DIR $MULTIQC_DIR

THREADS=8

echo "========================================"
echo " Step 1: Running FastQC on selected reads"
echo " Input directory: $SELECTED_DIR"
echo " Output directory: $QC_DIR"
echo " Threads: $THREADS"
echo "========================================"

fastqc -t $THREADS $SELECTED_DIR/*.fastq.gz -o $QC_DIR

echo ">>> FastQC completed."
echo ">>> Output files:"
ls -lh $QC_DIR | grep '.zip\|.html'

echo ""
echo "========================================"
echo " Step 2: Running MultiQC summary"
echo " Input: $QC_DIR"
echo " Output: $MULTIQC_DIR"
echo "========================================"

multiqc $QC_DIR -o $MULTIQC_DIR

echo ">>> MultiQC completed."
echo ">>> Report generated at: $MULTIQC_DIR/multiqc_report.html"
echo "========================================"
