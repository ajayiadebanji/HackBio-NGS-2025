#!/bin/bash
# Script: 05_fastp.sh
# Purpose: Perform adapter trimming and quality filtering with fastp
# Why:    Clean data improves assembly and reduces false positives in AMR/toxin detection.

SELECTED_DIR="selected"
TRIM_DIR="trimmed"
QC_DIR="qc/trimmed_fastqc"
FASTP_REPORTS="qc/fastp_reports"
MULTIQC_DIR="qc/multiqc_trimmed"

mkdir -p $TRIM_DIR $QC_DIR $FASTP_REPORTS $MULTIQC_DIR

THREADS=4

echo "========================================"
echo " Step 6: Running fastp on selected samples"
echo "========================================"

while read sample; do
  echo ">>> Processing sample: $sample"
  
  fastp \
    -i ${SELECTED_DIR}/${sample}_1.fastq.gz \
    -I ${SELECTED_DIR}/${sample}_2.fastq.gz \
    -o ${TRIM_DIR}/${sample}_1.trim.fastq.gz \
    -O ${TRIM_DIR}/${sample}_2.trim.fastq.gz \
    -w $THREADS \
    --detect_adapter_for_pe \
    --cut_front --cut_tail \
    --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 50 \
    --html ${FASTP_REPORTS}/${sample}_fastp.html \
    --json ${FASTP_REPORTS}/${sample}_fastp.json
done < selected_samples.txt

echo ">>> fastp trimming completed."

echo "========================================"
echo " Running MultiQC on fastp outputs"
echo "========================================"

multiqc $FASTP_REPORTS -o $MULTIQC_DIR

echo ">>> MultiQC summary generated at $MULTIQC_DIR/multiqc_report.html"
