#!/bin/bash -l
#SBATCH --mail-type=all

module load bioinfo-tools samtools

set -eux

in=$1
shift
outdir=$1
shift
samtools_flags=$1
shift
grep_command=$@

fnam=$(basename ${in/.bam/})

samtools view -h $samtools_flags $in | grep -F $grep_command | samtools view -b | samtools sort -o $outdir/$fnam.BWAfilter.bam
samtools index $outdir/$fnam.BWAfilter.bam

