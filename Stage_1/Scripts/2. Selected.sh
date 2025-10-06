#!/bin/bash
# 2. create selected

mkdir -p selected/raw
while read sample; do
  # adjust suffix patterns if your files differ (_1/_2 or _R1/_R2)
  if [[ -f "$PWD/${sample}_1.fastq.gz" && -f "$PWD/${sample}_2.fastq.gz" ]]; then
    ln -s "$PWD/${sample}_1.fastq.gz" selected/raw/${sample}_1.fastq.gz
    ln -s "$PWD/${sample}_2.fastq.gz" selected/raw/${sample}_2.fastq.gz
  else
    echo "Missing pair for $sample - check filename pattern"
  fi
done < selected_samples.txt
