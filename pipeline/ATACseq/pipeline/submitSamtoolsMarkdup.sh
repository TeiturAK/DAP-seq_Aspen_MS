#!/bin/bash

## be verbose and print
set -eux

# functions
source ../UPSCb-common/src/bash/functions.sh
module load bioinfo-tools samtools

## variables
proj=u2021028
mail=teitur.ahlgren.kalman@umu.se

in=/BWA
out=/BWA_markdup

if [ ! -d $out ]; then
  mkdir -p $out
fi

for f in $(find $in -name "*.bam"); do
	fnam=$(basename $f).samtools-markdup
	sbatch -A $proj -t 3:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
	-J $fnam -p node -n 20 runSamtoolsMarkdup.sh $f $out
done
