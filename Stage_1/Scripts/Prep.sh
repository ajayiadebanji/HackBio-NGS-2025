#!/usr#!/usr/bin/env bash
# create_and_run_pipeline.sh
# Creates project structure, scripts, and runs the full SA Polony WGS analysis pipeline.
# Usage:
#   chmod +x create_and_run_pipeline.sh
#   ./create_and_run_pipeline.sh [--downsample N] [--threads T] [--no-run]
#
# Options:
#   --downsample N   Keep randomly N samples (paired reads preserved) -- helpful if system resources limited
#   --threads T      Number of threads to use (default 4)
#   --no-run         Create files only; do not run the heavy pipeline steps
#
# NOTE: This script expects bioinformatics tools installed in PATH:
#   fastp, fastqc, multiqc, spades.py, quast, blastn, makeblastdb, abricate, mlst, Rscript, wget, curl
# If you need, use conda/mamba to install them (environment YAML is created below).
set -euo pipefail

# Default config
DOWNSAMPLE=0
THREADS=4
NORUN=false
WORKDIR="${PWD}/sa_polony_project"
DOWNLOAD_HELPER_URL="https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh"

# Parse CLI args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --downsample) DOWNSAMPLE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --no-run) NORUN=true; shift ;;
    --help) echo "Usage: $0 [--downsample N] [--threads T] [--no-run]"; exit 0 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

echo "[*] Workdir: $WORKDIR"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Create directories
echo "[*] Creating project directories..."
dirs=(
  scripts env data/raw data/clean results/fastqc results/assembly results/abricate results/blast results/mlst results/summary notebooks figures
)
for d in "${dirs[@]}"; do mkdir -p "$d"; done

# Write environment YAML for reproducibility
cat > env/environment.yml <<'YML'
name: sa-polony
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - fastp
  - fastqc
  - multiqc
  - spades
  - quast
  - blast
  - abricate
  - mlst
  - mashtree
  - r-base
  - r-tidyverse
  - wget
  - curl
  - seqkit
  - pigz
  - parallel
YML
echo "[*] Wrote env/environment.yml"
