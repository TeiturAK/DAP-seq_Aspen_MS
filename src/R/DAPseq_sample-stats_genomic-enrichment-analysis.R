#' ---
#' title: "DAP-seq sample stats and genomic enrichment analysis"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' This script looks at the general stats of TFs and open chromatin and 
#' what types of features are enriched with peaks taking the size of the peaks and the features into account.

#' # Libraries
suppressPackageStartupMessages({
  library(systemPipeR)
  library(ggplot2)
  library(ggridges)
  library(rlist)
  library(stringr)
  library(DT)
})

moduleload("bioinfo-tools BEDTools")

#' # Data
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered" 
DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.narrowPeak$", full.names = TRUE, all.files = FALSE, recursive = FALSE)
names(DAPseq.peak_files) <- substr(basename(DAPseq.peak_files), start = 1, stop = 11)
ATACseq.peaks.path <- "/mnt/picea/home/tkalman/ATAC-seq/aspen/macs2_p-0.05_18_May_2022/P15258_202_S4_L001_trimmomatic.sorted.dups-removed.BWAfilter_peaks.narrowPeak"

TF.sample_info <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/sample-info/P11505_P13101_P17252_merged-replicates_sampleinfo.tsv")

promoter.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/Potra02_promoters2kb_annotations.bed"

# three_prime_UTR.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_3primeUTR.sorted.merged.bed"
# five_prime_UTR.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_5primeUTR.sorted.merged.bed"
# CDS.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_CDS.sorted.merged.bed"
# intron.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intron.sorted.merged.bed"
# intergenic.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_intergenic_annotations.bed"
# promoter.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/Potra02_promoters2kb_annotations.bed"
# pure_intergenic.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_pure_intergenic_annotations.bed"
# LTR.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_LTR_annotations.neighbours-merged.bed"
# repeat.annotations <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/annotations-for-DAPseq-analysis/Potra02_repeat_annotations.neighbours-merged.LTR-subtract.bed"

#' Exclude the non-merged replicates from this analysis so that we do not amplify any results based on these samples
replicate_1 <- c("P13101_1005", "P17252_1089", "P17252_1189", "P17252_1245")
replicate_2 <- c("P13101_1115", "P17252_1090", "P17252_1190", "P17252_1246")
replicate_3 <- c("P13101_1396", "P17252_1091", "P17252_1191", "P17252_1247")
replicate_4 <- c("P13101_1518", "P17252_1092", "P17252_1192", "P17252_1248")
replicates <- c(replicate_1, replicate_2, replicate_3, replicate_4)

DAPseq.peak_files <- DAPseq.peak_files[!(names(DAPseq.peak_files) %in% replicates)]

#' # General TF peak stats
#' What is the number of TF samples?
print (length(DAPseq.peak_files))

#' How many of the TFs created peaks?
DAPseq.sample.peak_count <- do.call(rbind, lapply(DAPseq.peak_files, function (x) {
  TF.sample <- substr(basename(x), start = 1, stop = 11)
  peak.count <- as.numeric(system(paste("cat", x, "| wc -l"), intern = TRUE))
  data.frame(TF.sample, peak.count)
}))

print (length(DAPseq.sample.peak_count$peak.count[which(DAPseq.sample.peak_count$peak.count > 0)]))

#' What was the total peak count?
print (sum(DAPseq.sample.peak_count$peak.count))

#' Min, Max and median TF peak count? 
TF.peak_count.overview <- data.frame("min.count" = min(DAPseq.sample.peak_count$peak.count[which(DAPseq.sample.peak_count$peak.count > 0)]),
                                     "max.count" = max(DAPseq.sample.peak_count$peak.count),
                                     "median.count" = median(DAPseq.sample.peak_count$peak.count))

print (TF.peak_count.overview)

#' # General TF peak stats after ATAC intersect
DAPseq.ATAC_intersect.sample.peak_count <- do.call(rbind, lapply(DAPseq.peak_files, function (x) {
  TF.sample <- substr(basename(x), start = 1, stop = 11)
  peak.count <- as.numeric(system(paste("cat", x, "| bedtools intersect -wa -u -a - -b", ATACseq.peaks.path, "| wc -l"), intern = TRUE))
  data.frame(TF.sample, peak.count)
}))

#' How many TF samples had at least one peak intersecting open chromatin?
print (length(DAPseq.ATAC_intersect.sample.peak_count$peak.count[which(DAPseq.ATAC_intersect.sample.peak_count$peak.count > 0)]))

#' What was the total intersect count?
print (sum(DAPseq.ATAC_intersect.sample.peak_count$peak.count))

#' Min, Max and median TF peak count? 
TF.ATAC_intersect.peak_count.overview <- data.frame("min.count" = min(DAPseq.ATAC_intersect.sample.peak_count$peak.count[which(DAPseq.ATAC_intersect.sample.peak_count$peak.count > 0)]),
                                                    "max.count" = max(DAPseq.ATAC_intersect.sample.peak_count$peak.count),
                                                    "median.count" = median(DAPseq.ATAC_intersect.sample.peak_count$peak.count))

print (TF.ATAC_intersect.peak_count.overview)

#' # TF peak count in promoter after ATAC intersect
DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count <- do.call(rbind, list.clean(lapply(DAPseq.peak_files, function (x) {
  TF.sample <- substr(basename(x), start = 1, stop = 11)
  promoter.IDs <- unique(system(paste("bedtools intersect -wa -a", x, "-b", ATACseq.peaks.path, "| bedtools intersect -wa -a", promoter.annotations, "-b - | awk '{print $4}'"), intern = TRUE))
  if (length(promoter.IDs) > 0) {
    data.frame(TF.sample, promoter.IDs)
    }
  })))

#' How many promoters were enriched with at least one TFBS in open chromatin?
print (length(unique(DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count$promoter.IDs)))

#' Min, Max and median intersect count for TF samples intersecting open chromatin in promoter regions? 
DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count.split <- split(DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count, f = DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count$TF.sample)

DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count.split_count <- data.frame("promoter.count" = as.numeric(unlist(lapply(DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count.split, NROW))))

TF.ATAC_intersect.promoter_intersect.peak_count.overview <- data.frame("min.count" = min(DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count.split_count$promoter.count),
                                                                       "max.count" = max(DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count.split_count$promoter.count),
                                                                       "median.count" = median(DAPseq.ATAC_intersect.promoter_intersect.sample.peak_count.split_count$promoter.count))

print (TF.ATAC_intersect.promoter_intersect.peak_count.overview)

#' # Genome enrichment analysis
#' Genome annotations 
#' feature_list <- list(three_prime_UTR.annotations, 
#'                      five_prime_UTR.annotations, 
#'                      CDS.annotations, 
#'                      intron.annotations, 
#'                      intergenic.annotations,
#'                      promoter.annotations,
#'                      pure_intergenic.annotations,
#'                      LTR.annotations,
#'                      repeat.annotations)
#' 
#' names(feature_list) <- c("3primeUTR",
#'                          "5primeUTR",
#'                          "CDS",
#'                          "intron",
#'                          "intergenic",
#'                          "promoter",
#'                          "pure_intergenic",
#'                          "LTR",
#'                          "repeat")
#' 
#' #' What are the sizes of the genomic features?
#' feature_lengths.df <- data.frame(do.call(rbind, lapply(names(feature_list), function (tmp.name) {
#'   tmp <- feature_list[[tmp.name]]
#'   tmp.length <- as.numeric(system(paste("bedtools sort -i", tmp, "| bedtools merge -i - | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'"), intern = TRUE))
#'   data.frame("feature" = tmp.name,
#'              "total_length" = tmp.length)
#' })))
#' 
#' ggplot(feature_lengths.df, aes(x = feature, y = total_length, fill = feature)) +
#'   geom_col() +
#'   ggtitle("Genomic feature sizes") +
#'   theme_minimal()
#' 
#' #' DAP-seq - ATAC-seq intersect
#' DAPseq_ATACseq.peak_intersect.df <- data.frame(do.call(rbind, str_split(system(paste(paste0("cat ", DAPseq.dir_path, "/*.narrowPeak"), "| bedtools intersect -wa -a - -b", ATACseq.peaks.path, "| bedtools intersect -C -a - -b", 
#'                                                                              three_prime_UTR.annotations, 
#'                                                                              five_prime_UTR.annotations, 
#'                                                                              CDS.annotations, 
#'                                                                              intron.annotations, 
#'                                                                              intergenic.annotations,
#'                                                                              promoter.annotations,
#'                                                                              pure_intergenic.annotations,
#'                                                                              LTR.annotations,
#'                                                                              repeat.annotations,
#'                                                                              "-names 3primeUTR 5primeUTR CDS intron intergenic promoter pure_intergenic LTR repeat"),
#'                                                                        intern = TRUE), pattern = "\t")))
#' 
#' DAPseq_ATACseq.peak_intersect.df$X12 <- as.numeric(DAPseq_ATACseq.peak_intersect.df$X12)
#' DAPseq_ATACseq.peak_intersect.feature_split <- split(DAPseq_ATACseq.peak_intersect.df, f = DAPseq_ATACseq.peak_intersect.df$X11)
#' 
#' DAPseq_ATACseq.peak_intersect.feature_count.df <- data.frame("count" = do.call(rbind, lapply(DAPseq_ATACseq.peak_intersect.feature_split, function (x) {
#'   sum(x$X12)
#' })))
#' 
#' DAPseq_ATACseq.peak_intersect.feature_count.df$feature <- rownames(DAPseq_ATACseq.peak_intersect.feature_count.df)
#' 
#' ggplot(DAPseq_ATACseq.peak_intersect.feature_count.df, aes(x = feature, y = count, fill = feature)) +
#'   geom_col() +
#'   ggtitle("DAP-ATAC: Genomic feature intersect") +
#'   theme_minimal()
#' 
#' #' What were the intersect counts in pure numbers?
#' print (DAPseq_ATACseq.peak_intersect.feature_count.df)

#' #' What do our numbers look like compared to a random distribution?
#' contiglengths <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/Ptremula_annotations/contiglengths"
#' 
#' N <- 5
#' shuffled_peak.intersect.df <- data.frame(do.call(rbind, lapply(1:N, function (i) {
#'   tmp.shuffled_peak_intersect.df <- data.frame(do.call(rbind, str_split(system(paste(paste0("cat ", DAPseq.dir_path, "/*.narrowPeak"), "| bedtools shuffle -seed", i, "-i - -g", contiglengths, "| bedtools intersect -wa -a - -b", ATACseq.peaks.path, "| bedtools intersect -C -a - -b", 
#'                                                                                      three_prime_UTR.annotations, 
#'                                                                                      five_prime_UTR.annotations, 
#'                                                                                      CDS.annotations, 
#'                                                                                      intron.annotations, 
#'                                                                                      intergenic.annotations,
#'                                                                                      promoter.annotations,
#'                                                                                      pure_intergenic.annotations,
#'                                                                                      LTR.annotations,
#'                                                                                      repeat.annotations,
#'                                                                                      "-names 3primeUTR 5primeUTR CDS intron intergenic promoter pure_intergenic LTR repeat"),
#'                                                                                intern = TRUE), pattern = "\t")))
#'   tmp.shuffled_peak_intersect.df$X12 <- as.numeric(tmp.shuffled_peak_intersect.df$X12)
#'   tmp.shuffled_peak_intersect.feature_split <- split(tmp.shuffled_peak_intersect.df, f = tmp.shuffled_peak_intersect.df$X11)
#'   
#'   tmp.shuffled_peak_intersect.feature_count.df <- data.frame("count" = do.call(rbind, lapply(tmp.shuffled_peak_intersect.feature_split, function (x) {
#'     sum(x$X12)
#'   })))
#'   
#'   tmp.shuffled_peak_intersect.feature_count.df$feature <- rownames(tmp.shuffled_peak_intersect.feature_count.df)
#'   tmp.shuffled_peak_intersect.feature_count.df
#' })))
#' 
#' #' Need to fix the names before comparison
#' DAPseq_ATACseq.peak_intersect.feature_count.df$feature <- as.factor(c("three_prime_UTR.annotations","five_prime_UTR.annotations","CDS.annotations","intron.annotations","intergenic.annotations","promoter.annotations","pure_intergenic.annotations","LTR.annotations","repeat.annotations"))
#' 
#' ggplot(shuffled_peak.intersect.df, aes(x = count, y = feature, fill = feature)) + 
#'   geom_density_ridges(rel_min_height = 0.01, alpha = 0.75) +
#'   geom_segment(data = DAPseq_ATACseq.peak_intersect.feature_count.df, aes(x = count, xend = count, y = as.numeric(feature), yend = as.numeric(feature) + .9),
#'                color = "red") +
#'   scale_y_discrete(expand = c(0.01, 0)) +
#'   ggtitle("Shuffled peak feature intersect count distributions compared to observed intersect count (red line)") +
#'   xlab("Intersect (count)") +
#'   theme_minimal()
