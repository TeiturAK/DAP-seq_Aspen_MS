#' ---
#' title: "DAP-seq - ATAC-seq peak subsets for motifs and networks"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' DAP-seq - ATAC-seq peak subsets for motifs and networks

#' # Libraries
suppressPackageStartupMessages({
  library(systemPipeR)
  library(rlist)
  library(ggplot2)
  library(stringr)
  
})

moduleload("bioinfo-tools BEDTools")

#' # Annotation data
TF.sample_info <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/sample-info/P11505_P13101_P17252_merged-replicates_sampleinfo.tsv")
potri_potra.df <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/potra_potri_best_diamond.tsv",
                             header = FALSE, col.names = c("potra02", "potri"))
potra_artha.df <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/potra_artha_best_diamond.tsv", 
                             header = FALSE, col.names = c("potra", "artha"))

promoter.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/Potra02_promoters2kb_annotations.bed"


#' # Intersect DAP-seq peaks with promoters to create subset needed for motif analysis and write to file
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered"
DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.narrowPeak$", full.names = TRUE, all.files = FALSE, recursive = FALSE)
promoter_peaks.dir <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_promoter-intersected-peaks_only-neg-control-filtered"
# system(paste("mkdir", promoter_peaks.dir), intern = TRUE)
# DAPseq.peaks_intersecting_promoters <- lapply(DAPseq.peak_files, function (tmp.path) {
#   tmp.name <- basename(tmp.path)
#   system(paste("bedtools intersect -wa -u -a", tmp.path, "-b", promoter.annotations, ">", paste(promoter_peaks.dir, tmp.name, sep = "/")), intern = TRUE)
# })

#' Doing one extra for the MEME run where clusters have not been removed
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_only-neg-control-blacklist-filtered"
promoter_peaks.dir <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_promoter-intersected-peaks_only-neg-control-filtered"
# system(paste("mkdir", promoter_peaks.dir), intern = TRUE)

# DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.narrowPeak$", full.names = TRUE, all.files = FALSE, recursive = FALSE)
# DAPseq.peaks_intersecting_promoters <- lapply(DAPseq.peak_files, function (tmp.path) {
#   tmp.name <- basename(tmp.path)
#   system(paste("bedtools intersect -wa -u -a", tmp.path, "-b", promoter.annotations, ">", paste(promoter_peaks.dir, tmp.name, sep = "/")), intern = TRUE)
# })
# 
# DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.bed$", full.names = TRUE, all.files = FALSE, recursive = FALSE)
# DAPseq.peaks_intersecting_promoters <- lapply(DAPseq.peak_files, function (tmp.path) {
#   tmp.name <- basename(tmp.path)
#   system(paste("bedtools intersect -wa -u -a", tmp.path, "-b", promoter.annotations, ">", paste(promoter_peaks.dir, tmp.name, sep = "/")), intern = TRUE)
# })


#' # Intersect DAP-seq peaks with ATAC-seq peaks and then promoters to create global DAP+ATAC network inferred by promoter peaks
DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.narrowPeak$", full.names = TRUE, all.files = FALSE, recursive = FALSE)
ATACseq.peaks.path <- "/mnt/picea/home/tkalman/ATAC-seq/aspen/macs2_p-0.05_18_May_2022/P15258_202_S4_L001_trimmomatic.sorted.dups-removed.BWAfilter_peaks.narrowPeak"

DAPseq_ATACseq.peaks_intersecting_promoters.df <- data.frame(do.call(rbind, str_split(system(paste("cat", paste0(DAPseq.dir_path, "/*.narrowPeak"), "| bedtools sort -i - | bedtools intersect -wa -u -a - -b", ATACseq.peaks.path, "| bedtools intersect -wo -a - -b", promoter.annotations), intern = TRUE), pattern = "\t")))

DAPseq_ATACseq.peaks_intersecting_promoters.df$X4 <- substr(DAPseq_ATACseq.peaks_intersecting_promoters.df$X4, start = 1, stop = 11)
DAPseq_ATACseq.peaks_intersecting_promoters.df <- DAPseq_ATACseq.peaks_intersecting_promoters.df[, c("X1", "X2", "X3", "X4", "X14")]
colnames(DAPseq_ATACseq.peaks_intersecting_promoters.df) <- c("Seq", "Start", "End", "DAPseq.TF", "Potra02.Target")

DAPseq_ATACseq.peaks_intersecting_promoters.df$Potra02.TF <- TF.sample_info$Potra02[match(DAPseq_ATACseq.peaks_intersecting_promoters.df$DAPseq.TF, TF.sample_info$Sample.Name)]

DAPseq_ATACseq.peaks_intersecting_promoters.df <- DAPseq_ATACseq.peaks_intersecting_promoters.df[, c("Seq", "Start", "End", "DAPseq.TF", "Potra02.TF", "Potra02.Target")]

# write.table(DAPseq_ATACseq.peaks_intersecting_promoters.df, file = "/mnt/picea/home/tkalman/EvoTree/TF-target_networks/DAPseq_ATAC-filtered-peaks_intersect_w_promoter2kb_20-May-2022.tsv",
#              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
DAPseq_ATACseq.peaks_intersecting_promoters.df <- read.delim("/mnt/picea/home/tkalman/EvoTree/TF-target_networks/DAPseq_ATAC-filtered-peaks_intersect_w_promoter2kb_20-May-2022.tsv")

#' # Find TFs and targets are are of particular interest for analysis
TFs_of_interest.df <- openxlsx::read.xlsx("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/other-publications-for-validation/For Teitur 25.11.2021.xlsx")

TFs_of_interest.df$Potra02 <- potri_potra.df$potra02[match(TFs_of_interest.df$P..trichocarpa, potri_potra.df$potri)]
TFs_of_interest.df <- TFs_of_interest.df[!is.na(TFs_of_interest.df$Potra02), ]

TFs_of_interest.network.df <- DAPseq_ATACseq.peaks_intersecting_promoters.df[DAPseq_ATACseq.peaks_intersecting_promoters.df$Potra02.TF %in% TFs_of_interest.df$Potra02, ]

TFs_of_interest.network.df <- merge(TFs_of_interest.df, TFs_of_interest.network.df, by.x = "Potra02", by.y = "Potra02.TF")

TFs_of_interest.network.df <- TFs_of_interest.network.df[, c("DAPseq.TF", "Potra02", "P..trichocarpa", "Potri.name", "ATG", "Gene.name.in.Arabidopsis", "X7", "Potra02.Target")]

TFs_of_interest.network.df$Potri.Target <- potri_potra.df$potri[match(TFs_of_interest.network.df$Potra02.Target, potri_potra.df$potra02)]

TFs_of_interest.network.df$Arabidopsis.Target <- potra_artha.df$artha[match(TFs_of_interest.network.df$Potra02.Target, potra_artha.df$potra)]

colnames(TFs_of_interest.network.df) <- c("DAPseq.TF", "Potra02.TF", "P.trichocarpa.TF", "Comment: Potri.name", "Comment: ATG", "Comment: Gene.name.in.Arabidopsis", "Comment", "Potra02.Target", "Potri.Target", "Arabidopsis.Target")

# write.table(TFs_of_interest.network.df, file = "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/TFs_of_interest.network.tsv",
#              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
# openxlsx::write.xlsx(TFs_of_interest.network.df, file = "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/TFs_of_interest.network.xlsx")
