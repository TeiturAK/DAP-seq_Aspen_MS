# DAP-seq_Aspen_MS

## src/R
Scripts for analysis of DAP-seq and ATAC-seq data; creation of peak inferred regulatory network.

## src/py
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


Output is found in:
`/mnt/picea/home/tkalman/FOR-ANALYSIS-IN-ALL-PHD-PROJECTS/genomic_annotations/Ptremula_annotations/Potra02_genes.TSS.sorted.bed`

### Promoter annotations 
Created using `src/py/extract_promoter_V4.py` with gene annotations in `/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff` and 2000bp upstream as the length parameter.

Output is found in:
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/Potra02_promoters2kb_annotations.bed`

### Intergenic annotations 
Intergenic annotations were annotated by subtracting all intragenic regions from the contiglengts with the following steps.

Contiglengths were extracted from the first two fields of `/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_genome.fasta.fai`:
`awk '{print $1"\t"$2}' Potra02_genome.fasta.fai > contiglengths.bed`
Intragenic regions were extracted from the gene annotations in `/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff` and merged:
`grep gene Potra02_genes.gff | bedtools sort -i - | bedtools merge -i > intragenic.bed`

Lastly, intergenic regions were annotated by subtracting the intragenic from contiglengths with bedtools subtract:
`bedtools subtract -a contiglengths.bed -b intragenic.bed > intergenic`

Output is found in:
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intergenic_annotations.bed`

### LTR and repeat annotations 


Output is found in:
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_LTR_annotations.neighbours-merged.bed`
`/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_repeat_annotations.neighbours-merged.LTR-subtract.bed`
