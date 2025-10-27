# RNA-Seq analysis of _Culex pipiens_

- Author: Mia Scott
- Date: 10/24/25
- GitHub repository: exercise9
- This was run at the pitzer cluster of the Ohio Supercomputer Center (<www.osc.edu>)
- OSC working dir: /fs/ess/PAS2880/users/scott2291/week09/exercises

## Prerequisites

- **Input files**:
  - The input files are:
    - Paired-end FASTQ files
    - Reference genomes files (assembly and annotation)

- **Software**:
  - The scripts will run programs using Singularity/Apptainer containers
    available online, with links contained inside the scripts
  - If running at OSC, no software installations are needed
  - If you're not at OSC, you need to have `apptainer` and `Slurm` installed.

- You don't need to create the output dirs because the scripts will do so.

## Set up: Get the input files

1. Download the input FASTQ files and the reference genome annotation GTF
   from a GitHub repo as follows:

   ```bash
   git clone https://github.com/jelmerp/garrigos-data
   ```

2. Download the reference genome assembly FASTQ from NCBI as follows:

   ```bash
   URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/801/865/GCF_016801865.2_TS_CPP_V2/GCF_016801865.2_TS_CPP_V2_genomic.fna.gz
   wget -P garrigos-data/ref "$URL"
   ```

3. Unzip the reference genome files:

   ```bash
   gunzip -v garrigos-data/ref/*
   ```

## Set up: Set file paths

```bash
# Inputs:
fastq_dir=garrigos-data/fastq
ref_assembly=garrigos-data/ref/GCF_016801865.2_TS_CPP_V2_genomic.fna
ref_annotation=garrigos-data/ref/GCF_016801865.2.gtf

# Base output dir:
outdir=results
```

## Step 1: Read QC with FastQC

The `fastqc.sh` script runs on 1 FASTQ file at a time and takes 2 arguments -
the input FASTQ file and output dir:

```bash
for fastq in "$fastq_dir"/*fastq.gz; do
    sbatch scripts/fastqc.sh "$fastq" "$outdir"/fastqc
done
```

## Step 2: FastQC summaries with MultiQC

```bash
sbatch scripts/multiqc.sh "$outdir"/fastqc "$outdir"/multiqc
```

## Step 3: Read trimming with TrimGalore

```bash
for R1 in "$fastq_dir"/*_R1.fastq.gz; do
    R2="${R1/_R1/_R2}"
    sbatch scripts/trimgalore.sh "$R1" "$R2" "$outdir"/trimgalore
done

```

## Step 4: Creating a reference genome index with STAR

```bash
sbatch scripts/star_index.sh "$ref_assembly" "$ref_annotation" results/star/index
```

## Step 5: Trimmed read alignment to the reference genome index with STAR

```bash
for R1 in results/trimgalore/*R1.fastq.gz; do
    R2="${R1/_R1/_R2}"
    sbatch scripts/star_align.sh "$R1" "$R2" results/star/index "$ref_annotation" results/star
done
```