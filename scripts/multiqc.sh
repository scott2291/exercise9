#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-multiqc-%j.out
set -euo pipefail

# Constants
MULTIQC_CONTAINER=oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51

# Copy the placeholder variables
indir="$1"
outdir="$2"

# Initial logging
echo "# Starting script multiqc.sh"
date
echo "# Input dir:                      $indir"
echo "# Output file:                    $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run MultiQC
apptainer exec "$MULTIQC_CONTAINER" multiqc \
    --outdir "$outdir" \
    "$indir"

# Final logging
echo
echo "# Used MultiQC version:"
apptainer exec "$MULTIQC_CONTAINER" multiqc --version
echo "# Successfully finished script multiqc.sh"
date
