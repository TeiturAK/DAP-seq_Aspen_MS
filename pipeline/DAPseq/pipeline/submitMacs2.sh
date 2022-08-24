#!/bin/bash

## be verbose and print
set -eux

## source functions
source ../UPSCb-common/src/bash/functions.sh

## variables
proj=u2019003
mail=teitur.ahlgren.kalman@umu.se

in=/BWA_filtered

cont=/neg_control/neg_control.sorted.merged.bam

#P.tremula with mitochondria and chloroplast
effective_genome_size=382912866

#out=/macs2_cutoff-analysis
out=/macs2_p-0.001_15_May_2022

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## Execute
#for f in $(find $in -name "*.bam"); do
#  fnam=$(basename ${f/.bam/})
#  sbatch -A $proj -t 3:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
#  -J $fnam -p node -n 20 ../UPSCb-common/pipeline/runMacs2.sh -p -g $effective_genome_size -n $fnam $f $cont $out --cutoff-analysis
#done

for f in $(find $in -name "*.bam"); do
  fnam=$(basename ${f/.bam/})
  sbatch -A $proj -t 3:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
  -J $fnam -p node -n 20 ../UPSCb-common/pipeline/runMacs2.sh -p -g $effective_genome_size -n $fnam $f $cont $out -p 0.001
done

