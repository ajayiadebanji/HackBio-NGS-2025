nano bipolar_rnaseq.sh

#!/bin/bash
# ==============================================================
# 🧬 RNA-seq Pipeline: Unmasking Molecular Signatures of Bipolar II Disorder
# Author: Ajayi Adebanji
# Institution: HackBio Internship Program
# Date: October 2025
# ==============================================================
# Description:
# This pipeline performs RNA-seq analysis of iPSC-derived microglia 
# from Bipolar II disorder patients (familial and sporadic) and controls.
# Tools used: SRA Toolkit, FastQC, fastp, STAR, Samtools, featureCounts, R (DESeq2, clusterProfiler)
# ==============================================================

# ========== 1. Set up working environment ==========
WORKDIR=$HOME/BipolarII_RNAseq
mkdir -p $WORKDIR/{raw_data,trimmed,qc_reports,aligned,counts,genome,results}
cd $WORKDIR

# ========== 2. Download SRA data ==========
# Replace SRR IDs with actual dataset if needed
SRR_IDS=(
SRR33243164 SRR33243165 SRR33243166 SRR33243167
SRR33243168 SRR33243169 SRR33243170 SRR33243171
)

echo "📥 Downloading RNA-seq reads..."
for id in "${SRR_IDS[@]}"; do
  fasterq-dump $id -O raw_data --split-files
done

# ========== 3. Run FastQC on raw reads ==========
echo "🧪 Running FastQC on raw reads..."
mkdir -p qc_reports/raw
fastqc raw_data/*.fastq -o qc_reports/raw
multiqc qc_reports/raw -o qc_reports/raw

# ========== 4. Quality trimming with fastp ==========
echo "✂️ Trimming adapters and low-quality reads..."
mkdir -p trimmed qc_reports/trimmed
for id in "${SRR_IDS[@]}"; do
  fastp -i raw_data/${id}_1.fastq -I raw_data/${id}_2.fastq \
        -o trimmed/${id}_1.trimmed.fastq -O trimmed/${id}_2.trimmed.fastq \
        -h qc_reports/trimmed/${id}.html -j qc_reports/trimmed/${id}.json \
        -q 20 -u 30 -n 5 -l 30
done

fastqc trimmed/*.fastq -o qc_reports/trimmed
multiqc qc_reports/trimmed -o qc_reports/trimmed

# ========== 5. Download Human Reference Genome ==========
cd genome
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip *.gz

# ========== 6. Build STAR genome index ==========
mkdir -p STAR_index
STAR --runThreadN 8 --runMode genomeGenerate \
     --genomeDir STAR_index \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.112.gtf \
     --sjdbOverhang 100

cd $WORKDIR

# ========== 7. Align reads with STAR ==========
echo "🚀 Running STAR alignment..."
mkdir -p aligned
for id in "${SRR_IDS[@]}"; do
  STAR --runThreadN 8 \
       --genomeDir genome/STAR_index \
       --readFilesIn trimmed/${id}_1.trimmed.fastq trimmed/${id}_2.trimmed.fastq \
       --outFileNamePrefix aligned/${id}_ \
       --outSAMtype BAM SortedByCoordinate
done

# ========== 8. Index BAM files and get stats ==========
for bam in aligned/*.bam; do
  samtools index $bam
  samtools flagstat $bam > ${bam%.bam}.stats.txt
done

# ========== 9. Generate counts with featureCounts ==========
echo "📊 Quantifying gene expression..."
featureCounts -T 8 -p -B -C \
  -a genome/Homo_sapiens.GRCh38.112.gtf \
  -o counts/gene_counts.txt aligned/*.bam

# ========== 10. Run DESeq2 R script ==========
echo "🎯 Proceeding to DESeq2 differential analysis..."
Rscript DE_analysis_BipolarII.R

echo "✅ Pipeline complete! Results saved in /results"

# Making it executable
chmod +x bipolar_rnaseq.sh
./bipolar_rnaseq.sh


# Create the R script that performs DESeq2 analysis, visualization, and enrichment

nano DE_analysis_BipolarII.R
# ==============================================================
# 🎯 Differential Expression Analysis (DESeq2)
# Dataset: Bipolar II Disorder vs Control (iPSC-derived Microglia)
# Author: Ajayi Adebanji
# ==============================================================

# Load libraries
suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ========== 1. Import count data ==========
counts <- read.delim("counts/gene_counts.txt", comment.char="#", row.names=1)
counts <- counts[,6:ncol(counts)]  # Remove annotation columns
colnames(counts) <- c("Fam1","Fam2","Fam3","Fam4","Spo1","Spo2","Spo3","Spo4")

# Define sample metadata
condition <- factor(c("Familial","Familial","Familial","Familial",
                      "Sporadic","Sporadic","Sporadic","Sporadic"))
coldata <- data.frame(row.names=colnames(counts), condition)

# ========== 2. Create DESeq2 dataset ==========
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# ========== 3. Run DESeq2 ==========
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef="condition_Sporadic_vs_Familial", type="apeglm")

# Save results
write.csv(as.data.frame(res), "results/DESeq2_results.csv")

# ========== 4. PCA Plot ==========
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Bipolar II Disorder RNA-seq") +
  theme_minimal()
ggsave("results/PCA_plot.png")

# ========== 5. Volcano Plot ==========
EnhancedVolcano(res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'Differential Gene Expression: Sporadic vs Familial Bipolar II',
  pCutoff = 0.05, FCcutoff = 1.5,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2')
)
ggsave("results/Volcano_plot.png")

# ========== 6. Heatmap ==========
topGenes <- head(order(res$padj), 50)
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=coldata, show_rownames=FALSE,
         main="Top 50 Differentially Expressed Genes")
ggsave("results/Heatmap_top50.png")

# ========== 7. Functional Enrichment ==========
sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
ego <- enrichGO(gene=sig_genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05)
write.csv(as.data.frame(ego), "results/GO_Enrichment.csv")

# KEGG Pathway
ekegg <- enrichKEGG(gene=sig_genes, organism='hsa', pvalueCutoff=0.05)
write.csv(as.data.frame(ekegg), "results/KEGG_Enrichment.csv")

# Enrichment barplot
barplot(ego, showCategory=10, title="Top GO Biological Processes")
ggsave("results/GO_barplot.png")

# ========== 8. Save session ==========
save.image("results/DESeq2_session.RData")

cat("\n✅ DESeq2 analysis completed successfully!\nResults saved in /results\n")


