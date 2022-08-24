#!/bin/bash

## be verbose and print
set -eux

## source functions
source ../UPSCb-common/src/bash/functions.sh

## variables
proj=u2019003
mail=teitur.ahlgren.kalman@umu.se

in=/BWA_filtered
out=/macs2_p-0.05

#P.tremula with mitochondria and chloroplast
effective_genome_size=382912866

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute; WITHOUT a control sample
for f in $(find $in -name "*.bam"); do
  fnam=$(basename ${f/.bam/})
  sbatch -A $proj -t 24:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
  -J $fnam -p node ../UPSCb-common/pipeline/runMacs2.sh -c -p -g $effective_genome_size $f $out -p 0.05
done
