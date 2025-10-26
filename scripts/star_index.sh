#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-star_index-%j.out
set -euo pipefail

# Constants
STAR_CONTAINER=oras://community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4

# Copy the placeholder variables
fasta="$1"
gtf="$2"
outdir="$3"

# Initial logging
echo "# Starting script star_index.sh"
date
echo "# Input assembly FASTA file:      $fasta"
echo "# Input annotation GTF file:      $gtf"
echo "# Output dir:                     $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run STAR
apptainer exec "$STAR_CONTAINER" STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles "$fasta" \
    --sjdbGTFfile "$gtf" \
    --genomeDir "$outdir" \
    --runThreadN 16

#? - Note: This command was kept as simple as possible - if using STAR in your research,
#?   then look into the --sjdbOverhang and --genomeSAindexNbases options as well.
#? - Explanation of some of the options used:
#?   --runMode genomeGenerate  => Instead of aligning, this "run mode" generates a genome index

# Final logging
echo
echo "# Used STAR version:"
apptainer exec "$STAR_CONTAINER" STAR --version
echo "# Successfully finished script star_index.sh"
date
