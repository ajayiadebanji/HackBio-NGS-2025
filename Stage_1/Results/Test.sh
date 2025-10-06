#!/usr/bin/bash
#4. trim_qc.sh
# Description: Quality control on trimmed data

TRIMMED_DATA_DIR="./selected/trimmed"
QC_DIR="$TRIMMED_DATA_DIR/trimmed_qc"

# Create QC output directory
mkdir -p "$QC_DIR"
echo "Creating QC directory: $QC_DIR"

# Check if there are trimmed FASTQ files
if [ -z "$(ls -A $TRIMMED_DATA_DIR/*_*.trim.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed FASTQ files found in $TRIMMED_DATA_DIR"
    exit 1
fi

echo "=== Running FastQC on trimmed data ==="
fastqc -o "$QC_DIR" "$TRIMMED_DATA_DIR"/*_*.trim.fastq.gz

echo "=== Running MultiQC on trimmed QC reports ==="
multiqc "$QC_DIR" -o "$QC_DIR"

echo "QC reports saved to: $QC_DIR"
