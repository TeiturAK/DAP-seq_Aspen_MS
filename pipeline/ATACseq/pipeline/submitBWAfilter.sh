#!/bin/bash

## be verbose and print
set -eux

# functions
source ../UPSCb-common/src/bash/functions.sh

## variables
proj=snic2021-22-932
mail=teitur.ahlgren.kalman@umu.se

in=/BWA_markdup
out=/BWA_filtered

# aspen
grep_command='grep -v -e KP861984.1 -e KT337313.1'


if [ ! -d $out ]; then
  mkdir -p $out
fi

for f in $(find $in -name "*.dups-removed.bam"); do
        fnam=$(basename $f)
        sbatch -A $proj -t 3:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
        -J $fnam -p node runBWAfilter.sh $f $out $grep_command
done
