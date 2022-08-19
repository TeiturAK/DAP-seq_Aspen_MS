# DAP-seq_Aspen_MS

## Content of src/R
Scripts for analysis of DAP-seq and ATAC-seq data; creation of peak inferred regulatory network.

## Content of src/py
Script for extracing promoter region from gene annotation file.

## Notes on the genomic feature annotations
### Intragenic annotations
5'UTR, 3'UTR, intron and CDS annotations where all extracted from the gene annotation file in `/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff` using the following grep and bedtools command:
`grep <genomic feature> Potra02_genes.gff | bedtools sort -i - | bedtools merge -s -i - > <genomic feature sorted and overlaps merged>` 

The created output is found in:
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_3primeUTR.sorted.merged.bed`
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_5primeUTR.sorted.merged.bed`
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_CDS.sorted.merged.bed`
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intron.sorted.merged.bed`

### TSS annotation
TSS annotations were created by extracting the first coordinate of every gene in the gene annotation file `/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff` with the following commands:
`grep gene Potra02_genes.gff | awk '{if($7== "+") print $1"\t"$4"\t"$4+1}' > TSS.plus-strand.bed`
`grep gene Potra02_genes.gff | awk '{if($7== "-") print $1"\t"$5-1"\t"$5}' > TSS.minus-strand.bed`
`cat TSS.plus-strand.bed TSS.minus-strand.bed | bedtools sort -i - > Potra02_genes.TSS.sorted.bed`

Output is found in:
`/mnt/picea/home/tkalman/FOR-ANALYSIS-IN-ALL-PHD-PROJECTS/genomic_annotations/Ptremula_annotations/Potra02_genes.TSS.sorted.bed`

### Promoter annotations 
Created using `src/py/extract_promoter_V4.py` with gene annotations in `/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff` and `2000` given as the upstream cutoff length parameter.

Output is found in:
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/Potra02_promoters2kb_annotations.bed`

### Intergenic annotations 
Intergenic annotations were annotated by subtracting all intragenic regions from the contiglengts with the following steps.

Contiglengths were extracted from the first two fields of `/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_genome.fasta.fai`:
`awk '{print $1"\t"$2}' Potra02_genome.fasta.fai > contiglengths.bed`
Intragenic regions were extracted from the gene annotations in `/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff` and merged:
`grep gene Potra02_genes.gff | bedtools sort -i - | bedtools merge -i > intragenic.bed`

Lastly, intergenic regions were annotated by subtracting the intragenic from contiglengths with bedtools subtract:
`bedtools subtract -a contiglengths.bed -b intragenic.bed > Potra02_intergenic_annotations.bed`

Output is found in:
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intergenic_annotations.bed`

### LTR and repeat annotations 



Repeats not overlapping annotated LTR were extracted from `/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_repeatmasked.gff.gz` with the following commands:
`zcat Potra02_repeatmasked.gff.gz | bedtools sort -i - | bedtools merge -s -i - > Potra02_repeat_annotations.neighbours-merged.bed`
`bedtools subtract -a Potra02_repeat_annotations.neighbours-merged.bed -b Potra02_LTR_annotations.neighbours-merged.bed > Potra02_repeat_annotations.neighbours-merged.LTR-subtract.bed`

Output is found in:
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_LTR_annotations.neighbours-merged.bed`
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_repeat_annotations.neighbours-merged.LTR-subtract.bed`

#### Some additional notes on the original files as there was uncertainty at one point

Repeat source

LTR source
LTR annotation was created by Bastian with [GitHub - oushujun/LTR_retriever: LTR_retriever is a highly accurate and sensitive program for identification of LTR retrotransposons; The LTR Assembly Index (LAI) is also included in this package.](https://github.com/oushujun/LTR_retriever#inputs). 

