#!/bin/bash

## be verbose and print
set -eux

## variables
proj=u2019003
mail=teitur.kalman@umu.se

in=/BWA
out=/BWA_filtered

if [ ! -d $out ]; then
  mkdir -p $out
fi

samtools_flags="-f 3 -F 12 -F 512 -q 30"
grep_command="-v -e KP861984.1 -e KT337313.1"

for f in $(find $in -name "*.bam"); do
	fnam=$(basename ${f/.bam/})
	sbatch -A $proj -t 3:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
	-J $fnam -p node -n 20 runBWAfilter.sh $f $out "$samtools_flags" "$grep_command"
done
