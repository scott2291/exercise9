#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-fastqc-%j.out
set -euo pipefail

# Load the OSC module for FastQC
module load fastqc/0.12.1

# Copy the placeholder variables
fastq="$1"
outdir="$2"

# Initial logging
echo "# Starting script fastqc.sh"
date
echo "# Input FASTQ file:   $fastq"
echo "# Output dir:         $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run FastQC
fastqc \
    --threads 4 \
    --outdir "$outdir" \
    "$fastq"

# Final logging
echo
echo "# Used FastQC version:"
fastqc --version
echo
echo "# Successfully finished script fastqc.sh"
date
