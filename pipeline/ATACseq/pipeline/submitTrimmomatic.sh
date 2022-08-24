#!/bin/bash

## be verbose and print
set -eux

export SINGULARITY_BINDPATH="/mnt:/mnt"

proj=u2021028
mail=teitur.kalman@umu.se

in=/raw
out=/trimmomatic

adpt=/mnt/picea/Modules/apps/bioinfo/trimmomatic/0.39/adapters/NexteraPE-PE.fa
clip="ILLUMINACLIP:$adpt:2:30:10:1:TRUE"
trim="SLIDINGWINDOW:5:20 MINLEN:50"
simg=/mnt/picea/projects/singularity/kogia/trimmomatic_0.39.sif

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute
for f in $(find $in -type f -name "*_R1_001.fastq.gz"); do
  fnam=$(basename ${f/.fastq.gz/})
  sbatch -A $proj --mail-user=$mail -t 01:00:00 -e $out/$fnam.err -o $out/$fnam.out \
  -J $fnam -p node ../UPSCb-common/pipeline/runTrimmomatic.sh -c $clip $simg $adpt $f ${f/_R1_001.fastq.gz/_R2_001.fastq.gz} $out $trim
done



