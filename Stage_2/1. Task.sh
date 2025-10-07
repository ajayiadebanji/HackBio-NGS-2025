#!/bin/bash
###############################################################################
# Project: Unmasking Molecular Signatures of Bipolar II Disorder
# Author: Ajayi Adebanji
# Date: October 2025
# Pipeline: RNA-seq analysis (STAR + featureCounts + DESeq2)
#
# Description:
# This pipeline performs end-to-end RNA-seq preprocessing and quantification:
# 1. Download FASTQ data from SRA
# 2. Perform quality control (FastQC + MultiQC)
# 3. Trim adapters and low-quality bases (Trim Galore)
# 4. Align reads to the human GRCh38 reference genome using STAR
# 5. Generate gene-level count matrices using featureCounts
# 6. Prepare for DESeq2 differential expression analysis in R
#
# Software required:
#   - SRA Toolkit
#   - FastQC
#   - MultiQC
#   - Trim Galore
#   - STAR
#   - Subread (featureCounts)
#   - R with DESeq2, pheatmap, EnhancedVolcano, clusterProfiler
###############################################################################

# ===========================
# 0. Set up directories
# ===========================
mkdir -p data raw trimmed fastqc multiqc star_index alignments counts logs

# ===========================
# 1. Download raw FASTQ files
# ===========================
echo "=== Downloading RNA-seq data from SRA ==="
cd data
for srr in SRR33243164 SRR33243165 SRR33243166 SRR33243167 SRR33243168 SRR33243169 SRR33243170 SRR33243171
do
    fasterq-dump $srr -O ../raw --split-files --threads 6
done
cd ..

# ===========================
# 2. Quality control (raw)
# ===========================
echo "=== Running FastQC on raw reads ==="
fastqc raw/*.fastq -o fastqc
multiqc fastqc -o multiqc/raw_qc_report

# ===========================
# 3. Adapter trimming
# ===========================
echo "=== Trimming adapters and low-quality reads ==="
for file in raw/*_1.fastq
do
    base=$(basename $file _1.fastq)
    trim_galore --paired raw/${base}_1.fastq raw/${base}_2.fastq -o trimmed
done

# QC after trimming
fastqc trimmed/* -o fastqc
multiqc fastqc -o multiqc/trimmed_qc_report

# ===========================
# 4. Genome indexing and alignment
# ===========================
echo "=== Building STAR index and aligning reads ==="

# Download GRCh38 reference if not already present
if [ ! -f "star_index/Genome" ]; then
    wget -O GRCh38.fa.gz https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip GRCh38.fa.gz

    wget -O genes.gtf.gz https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
    gunzip genes.gtf.gz

    STAR --runThreadN 8 \
         --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles GRCh38.fa \
         --sjdbGTFfile genes.gtf \
         --sjdbOverhang 100
fi

# Alignment
for sample in trimmed/*_1_val_1.fq
do
    base=$(basename $sample _1_val_1.fq)
    STAR --genomeDir star_index \
         --readFilesIn trimmed/${base}_1_val_1.fq trimmed/${base}_2_val_2.fq \
         --runThreadN 8 \
         --outFileNamePrefix alignments/${base}_ \
         --outSAMtype BAM SortedByCoordinate
done

# ===========================
# 5. Counting features
# ===========================
echo "=== Counting aligned reads per gene ==="
featureCounts -T 8 -p -a genes.gtf -o counts/gene_counts.txt alignments/*_Aligned.sortedByCoord.out.bam

# ===========================
# 6. Organize output
# ===========================
mkdir -p results
cp counts/gene_counts.txt results/
echo "Pipeline complete. Proceed to R for DESeq2 analysis."
###############################################################################
