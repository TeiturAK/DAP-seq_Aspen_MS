# DAP-seq_Aspen_MS

## src/R
Scripts for analysis of DAP-seq and ATAC-seq data; creation of peak inferred regulatory network.

## src/py
Script for extracing promoter region from gene annotation file.

## Notes on the genomic feature annotations
5'UTR, 3'UTR, intron and CDS annotations where all extracted from the gene annotation file 
### three_prime_UTR.annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_3primeUTR.sorted.merged.bed"

### five_prime_UTR.annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_5primeUTR.sorted.merged.bed"

### CDS.annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_CDS.sorted.merged.bed"

### intron.annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intron.sorted.merged.bed"



### intergenic.annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intergenic_annotations.bed"

### promoter.annotations 
Created using src/py/extract_promoter_V4.py with gene annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/Potra02_promoters2kb_annotations.bed"
src/py

### LTR.annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_LTR_annotations.neighbours-merged.bed"

### repeat.annotations 
"/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_repeat_annotations.neighbours-merged.LTR-subtract.bed"
