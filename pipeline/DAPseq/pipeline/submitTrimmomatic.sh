#!/bin/bash

## be verbose and print
set -ex

## variables
proj=u2019003
mail=teitur.kalman@umu.se

## tools
module load bioinfo-tools trimmomatic

in=/raw
out=/trimmomatic

clip="ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE-2.fa:2:30:10:2:TRUE"
trim="SLIDINGWINDOW:5:20 MINLEN:35"

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## Execute
for f in $(find $in -type f -name "*_R1_001.fastq.gz"); do
  fnam=$(basename ${f/.fastq.gz/})
  sbatch -A $proj --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
  -J $fnam -p rbx ../UPSCb-common/pipeline/runTrimmomatic.sh -c $clip $f ${f/_R1_001.fastq.gz/_R2_001.fastq.gz} $out $trim
done

