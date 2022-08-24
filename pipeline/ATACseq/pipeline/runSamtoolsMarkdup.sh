#!/bin/bash -l
#SBATCH --mail-type=all

module load bioinfo-tools samtools

set -ex

infile=$1
outdir=$2

filename=$(basename ${infile/.sorted.bam/})

## from https://www.biostars.org/p/415831/

samtools sort -n -o $outdir/$filename.name-sorted.bam -O BAM $infile

samtools fixmate -m $outdir/$filename.name-sorted.bam $outdir/$filename.fixmate.bam

samtools sort -o $outdir/$filename.sorted.bam $outdir/$filename.fixmate.bam

samtools markdup -S -r -s $outdir/$filename.sorted.bam $outdir/$filename.sorted.dups-removed.bam

samtools index $outdir/$filename.sorted.dups-removed.bam

rm $outdir/$filename.name-sorted.bam
rm $outdir/$filename.fixmate.bam
rm $outdir/$filename.sorted.bam
