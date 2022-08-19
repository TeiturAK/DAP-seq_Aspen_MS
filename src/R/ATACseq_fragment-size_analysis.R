#' ---
#' title: "ATAC-seq fragment size analysis"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' Analyse ATAC-seq insert sizes

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/home/tkalman")
#' ```

suppressPackageStartupMessages({
  library(ATACseqQC)
})

#' # Aspen wood for DAP-seq project
ATACseqQC::fragSizeDist("/mnt/picea/home/tkalman/ATAC-seq/aspen/bwa_markdup_hardmasked-w-mitochondrion-and-chloroplast_only-unique-reads_1-Mar-2022/P15258_202_S4_L001_trimmomatic.sorted.unique_mapped.mt_cp_filtered.sorted.dups-removed.bam", 
                        bamFiles.labels = "Aspen Wood")