#!/bin/bash

## be verbose and print
set -eux

## source functions
source ../UPSCb-common/src/bash/functions.sh

## variables
proj=u2021028
mail=teitur.ahlgren.kalman@umu.se

in=/trimmomatic
out=/BWA
inx=/Potra02_genome_with_chloroplast_and_mitochondrion.fasta.gz

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi


## execute
for f in $(find $in -name "*_trimmomatic_1.fq.gz"); do
  fnam=$(basename ${f/_trimmomatic_1.fq/})
  sbatch -A $proj -t 03:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
  -J $fnam -p node -n 20 --mem=40GB ../UPSCb-common/pipeline/runBwamem.sh -t 20 -f $f -r ${f/_trimmomatic_1.fq.gz/_trimmomatic_2.fq.gz} -i $inx -o $out
done
