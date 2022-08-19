#' ---
#' title: "DAP-seq - genomic feature enrichment"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' Looking at what genomic features the DAP-seq peaks intersect

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/home/tkalman")
#' ```

#' # Libraries
suppressPackageStartupMessages({
  library(systemPipeR)
  library(ggplot2)
  library(stringr)
})

moduleload("bioinfo-tools BEDTools")

#' Peak data
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered" 
ATACseq.peaks.path <- "/mnt/picea/home/tkalman/ATAC-seq/aspen/macs2_p-0.05_18_May_2022/P15258_202_S4_L001_trimmomatic.sorted.dups-removed.BWAfilter_peaks.narrowPeak"

#' Detailed annotations 
three_prime_UTR.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_3primeUTR.sorted.merged.bed"
five_prime_UTR.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_5primeUTR.sorted.merged.bed"
CDS.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_CDS.sorted.merged.bed"
intron.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intron.sorted.merged.bed"
intergenic.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intergenic_annotations.bed"
promoter.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/Potra02_promoters2kb_annotations.bed"
LTR.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_LTR_annotations.neighbours-merged.bed"
repeat.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_repeat_annotations.neighbours-merged.LTR-subtract.bed"

#' Detail annotations intersect, scaled to feature length
feature_list <- list(three_prime_UTR.annotations,
                     five_prime_UTR.annotations,
                     CDS.annotations,
                     intron.annotations,
                     intergenic.annotations,
                     promoter.annotations,
                     LTR.annotations,
                     repeat.annotations)

names(feature_list) <- c("3primeUTR",
                         "5primeUTR",
                         "CDS",
                         "intron",
                         "intergenic",
                         "promoter",
                         "LTR",
                         "repeat")


feature_lengths.df <- data.frame(do.call(rbind, lapply(names(feature_list), function (tmp.name) {
  tmp <- feature_list[[tmp.name]]
  tmp.length <- as.numeric(system(paste("bedtools sort -i", tmp, "| bedtools merge -i - | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'"), intern = TRUE))
  data.frame("feature" = tmp.name,
             "total_length" = tmp.length)
})))

ggplot(feature_lengths.df, aes(x = feature, y = total_length, fill = feature)) +
  geom_col() +
  ggtitle("Genomic feature sizes") +
  theme_minimal()


#' # DAP-seq - ATAC-seq intersect
DAPseq_ATACseq.peak_intersect.df <- data.frame(do.call(rbind, str_split(system(paste(paste0("cat ", DAPseq.dir_path, "/*.narrowPeak"), "| bedtools intersect -wa -a - -b", ATACseq.peaks.path, "| bedtools intersect -C -a - -b",
                                                                             three_prime_UTR.annotations,
                                                                             five_prime_UTR.annotations,
                                                                             CDS.annotations,
                                                                             intron.annotations,
                                                                             intergenic.annotations,
                                                                             promoter.annotations,
                                                                             LTR.annotations,
                                                                             repeat.annotations,
                                                                             "-names 3primeUTR 5primeUTR CDS intron intergenic promoter LTR repeat"),
                                                                       intern = TRUE), pattern = "\t")))

DAPseq_ATACseq.peak_intersect.df$X12 <- as.numeric(DAPseq_ATACseq.peak_intersect.df$X12)
DAPseq_ATACseq.peak_intersect.feature_split <- split(DAPseq_ATACseq.peak_intersect.df, f = DAPseq_ATACseq.peak_intersect.df$X11)

DAPseq_ATACseq.peak_intersect.feature_count.df <- data.frame("count" = do.call(rbind, lapply(DAPseq_ATACseq.peak_intersect.feature_split, function (x) {
  sum(x$X12)
})))

DAPseq_ATACseq.peak_intersect.feature_count.df$feature <- as.factor(rownames(DAPseq_ATACseq.peak_intersect.feature_count.df))

DAPseq_ATACseq.peak_intersect.feature_count.df$feature <- factor(DAPseq_ATACseq.peak_intersect.feature_count.df$feature, 
                                                                 levels = c("intergenic", "promoter", "5primeUTR", "intron", "CDS", "3primeUTR", "LTR", "repeat"))


ggplot(DAPseq_ATACseq.peak_intersect.feature_count.df, aes(x = feature, y = count, fill = feature)) +
  geom_col() +
  ggtitle("DAP-ATAC: Genomic feature intersect") +
  theme_minimal()

#' What were the intersect counts in pure numbers?
print (DAPseq_ATACseq.peak_intersect.feature_count.df)

#' # Now taking the size of the features into account as well
feature_lengths.df <- data.frame(do.call(rbind, lapply(names(feature_list), function (tmp.name) {
  tmp <- feature_list[[tmp.name]]
  tmp.length <- as.numeric(system(paste("bedtools sort -i", tmp, "| bedtools merge -i - | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'"), intern = TRUE))
  data.frame("feature" = tmp.name,
             "total_length" = tmp.length)
})))

DAPseq_ATACseq.peak_intersect.feature_count.df <- merge(DAPseq_ATACseq.peak_intersect.feature_count.df, feature_lengths.df)
DAPseq_ATACseq.peak_intersect.feature_count.df$count_scaled_to_length <- DAPseq_ATACseq.peak_intersect.feature_count.df$count/DAPseq_ATACseq.peak_intersect.feature_count.df$total_length

ggplot(DAPseq_ATACseq.peak_intersect.feature_count.df, aes(x = feature, y = count_scaled_to_length, fill = feature)) +
  geom_col() +
  ggtitle("Feature count/Feature length") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


