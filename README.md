# DAP-seq_Aspen_MS

## Contents of pipeline
SBATCH scripts for processing raw ATAC/DAP-seq data with Trimmomatic, BWA-MEM, Samtools Markdup, filtering reads aligned to mitochondria/chloroplast sequences and calling peaks with MACS2.

## Contents of src/R
Scripts for analysis of DAP-seq and ATAC-seq data; creation of peak inferred regulatory network.

## Contents of src/py
Script for extracing promoter region from gene annotation file.


## Notes on the genomic feature annotations used in the analysis
Gene and repeat/LTR  annotations were all collected from: Schiffthaler B, Delhomme N, Bernhardsson C, Jenkins J, Jansson S, Ingvarsson P, Schmutz J, Street N. 2019. An improved genome assembly of the European aspen Populus tremula. bioRxiv: 805614.

### Intragenic annotations
5'UTR, 3'UTR, intron and CDS annotations where all extracted from the gene annotation `Potra02_genes.gff` using the following grep and bedtools command:
`grep <genomic feature> Potra02_genes.gff | bedtools sort -i - | bedtools merge -s -i - > <genomic feature sorted and overlaps merged>` 

### TSS annotation
TSS annotations were created by extracting the first coordinate of every gene in the gene annotation file `Potra02_genes.gff` with the following commands:
`grep gene Potra02_genes.gff | awk '{if($7== "+") print $1"\t"$4"\t"$4+1}' > TSS.plus-strand.bed`
`grep gene Potra02_genes.gff | awk '{if($7== "-") print $1"\t"$5-1"\t"$5}' > TSS.minus-strand.bed`
`cat TSS.plus-strand.bed TSS.minus-strand.bed | bedtools sort -i - > Potra02_genes.TSS.sorted.bed`

### Promoter annotations 
Created using `src/py/extract_promoter_V4.py` with gene annotations in `Potra02_genes.gff` and `2000` given as the upstream cutoff length parameter.

### Intergenic annotations 
Intergenic annotations were annotated by subtracting all intragenic regions from the contiglengts with the following steps.

Contiglengths were extracted from the first two fields of `Potra02_genome.fasta.fai`:
`awk '{print $1"\t"$2}' Potra02_genome.fasta.fai > contiglengths.bed`
Intragenic regions were extracted from the gene annotations in `Potra02_genes.gff` and merged:
`grep gene Potra02_genes.gff | bedtools sort -i - | bedtools merge -i > intragenic.bed`

Lastly, intergenic regions were annotated by subtracting the intragenic from contiglengths with bedtools subtract:
`bedtools subtract -a contiglengths.bed -b intragenic.bed > Potra02_intergenic_annotations.bed`

### LTR and repeat annotations 

LTR annotations were extracted from `Potra02_LTR.gff3.gz` using the following command:
`zcat Potra02_LTR.gff3.gz | bedtools sort -i - | bedtools merge -s -d 1 -i - > Potra02_LTR_annotations.neighbours-merged.bed`

Repeats not overlapping annotated LTR were extracted from `Potra02_repeatmasked.gff.gz` with the following commands:
`zcat Potra02_repeatmasked.gff.gz | bedtools sort -i - | bedtools merge -s -d 1 -i - > Potra02_repeat_annotations.neighbours-merged.bed`
`bedtools subtract -a Potra02_repeat_annotations.neighbours-merged.bed -b Potra02_LTR_annotations.neighbours-merged.bed > Potra02_repeat_annotations.neighbours-merged.LTR-subtract.bed`
