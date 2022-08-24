#!/bin/bash

## be verbose and print
set -eux

## variables
proj=u2019003
mail=teitur.kalman@umu.se

## tools
module load bioinfo-tools samtools

dir=/BWA_filtered

samples="\
$dir/P17252_1001_S1_L001_trimmomatic.sorted.BWAfilter.bam \
$dir/P17252_1188_S188_L002_trimmomatic.sorted.BWAfilter.bam \
$dir/P17252_1244_S148_L002_trimmomatic.sorted.BWAfilter.bam \
$dir/P13101_1087_S87_L001_trimmomatic.sorted.BWAfilter.bam \
$dir/P13101_1187_S187_L002_trimmomatic.sorted.BWAfilter.bam \
$dir/P13101_1471_S93_L001_trimmomatic.sorted.BWAfilter.bam \
$dir/P13101_1576_S96_L001_trimmomatic.sorted.BWAfilter.bam"

out=/neg_control/neg_control.sorted.merged.bam

outdir=$(dirname $out)

if [ ! -d $outdir ]; then
  mkdir -p $outdir
fi

fnam=$(basename ${out/.bam/})

#echo $samples
sbatch -A $proj -t 3:00:00 --mail-user=$mail -e $outdir/$fnam.err -o $outdir/$fnam.out \
-J $fnam -p core -n 20 ../UPSCb-common/pipeline/runSamtoolsMerge.sh $out $samples

