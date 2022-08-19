#' ---
#' title: "DAP-seq peak QC"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' For a given set of TFs with peaks: Check if there are any odd peak clusters.
#' Characterize suspicious peaks like overrepresentation in scaffolds.
#' Generate a bed filter file with regions that are not desirable to be kept in the analysis.

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/home/tkalman")
#' ```

#' # Libraries
suppressPackageStartupMessages({
  library(systemPipeR)
  library(rlist)
  library(ggplot2)
  library(stringr)
  
})

moduleload("bioinfo-tools BEDTools")

#' # Sample information
TF.sample_info <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/sample-info/P11505_P13101_P17252_merged-replicates_sampleinfo.tsv")

#' # Peak data
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_p-0.001_15_May_2022" 
DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.narrowPeak$", full.names = TRUE, all.files = FALSE, recursive = FALSE)
names(DAPseq.peak_files) <- substr(basename(DAPseq.peak_files), start = 1, stop = 11)

#' What is the peak count in the different samples?
DAPseq.sample.peak_count <- do.call(rbind, lapply(DAPseq.peak_files, function (x) {
  TF.sample <- substr(basename(x), start = 1, stop = 11)
  peak.count <- as.numeric(system(paste("cat", x, "| wc -l"), intern = TRUE))
  data.frame(TF.sample, peak.count)
}))

DAPseq.sample.peak_count$Run <- substr(DAPseq.sample.peak_count$TF.sample, start = 1, stop = 6)
DAPseq.sample.peak_count$Pool <- TF.sample_info$UDF.Pooling[match(DAPseq.sample.peak_count$TF.sample, TF.sample_info$Sample.Name)]

#' Combining run and pool to identify plates prepared together
DAPseq.sample.peak_count$Batch <- paste0(DAPseq.sample.peak_count$Run, "-", DAPseq.sample.peak_count$Pool)

#' Plotting peak count
#' Removing samples with zero observations 
DAPseq.sample.peak_count <- DAPseq.sample.peak_count[which(DAPseq.sample.peak_count$peak.count > 0), ]

#' How many samples had peaks?
print (NROW(DAPseq.sample.peak_count))

ggplot(DAPseq.sample.peak_count, aes(x = peak.count)) +
  geom_density() +
  ggtitle("Per sample peak count") +
  xlab("Peak (count)") +
  theme_minimal()

#' # Identifying clusters
#' Exclude the non-merged replicates from this analysis so that we do not identify any clusters due from intersecting replicate samples
replicate_1 <- c("P13101_1005", "P17252_1089", "P17252_1189", "P17252_1245")
replicate_2 <- c("P13101_1115", "P17252_1090", "P17252_1190", "P17252_1246")
replicate_3 <- c("P13101_1396", "P17252_1091", "P17252_1191", "P17252_1247")
replicate_4 <- c("P13101_1518", "P17252_1092", "P17252_1192", "P17252_1248")
replicates <- c(replicate_1, replicate_2, replicate_3, replicate_4)

DAPseq.peak_files <- DAPseq.peak_files[!(names(DAPseq.peak_files) %in% replicates)]

#' Peak cluster intersect
# system(paste("mkdir /mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis"), intern = TRUE)

# DAPseq.merged_peaks.df <- do.call(rbind, str_split(system(paste("cat", paste(DAPseq.peak_files, collapse = " "), "| bedtools sort -i - | bedtools merge -i -"), intern = TRUE), pattern = "\t"))
# write.table(DAPseq.merged_peaks.df, file = "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.all_peaks.sorted.merged.bed",
#              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
DAPseq.merged_peaks.path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.all_peaks.sorted.merged.bed"

# DAPseq.all_peaks.df <- do.call(rbind, str_split(system(paste("cat", paste(DAPseq.peak_files, collapse = " "), "| bedtools sort -i -"), intern = TRUE), pattern = "\t"))
# write.table(DAPseq.all_peaks.df, file = "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.all_peaks.sorted.bed",
#              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
DAPseq.all_peaks.path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.all_peaks.sorted.bed"


DAPseq.peak_intersect.df <- data.frame(do.call(rbind, str_split(system(paste("bedtools intersect -c -a", DAPseq.merged_peaks.path, "-b", DAPseq.all_peaks.path),
                                                     intern = TRUE), pattern = "\t")))

DAPseq.peak_intersect.df$X4 <- as.numeric(DAPseq.peak_intersect.df$X4)

ggplot(DAPseq.peak_intersect.df, aes(x = X4)) +
  geom_density() +
  ggtitle("Peaks intersecting") +
  xlab("count of peaks intersecting") +
  theme_minimal()

DAPseq.peak_intersect.below_thresh.df <- DAPseq.peak_intersect.df[which(DAPseq.peak_intersect.df$X4 < 100), ]

ggplot(DAPseq.peak_intersect.below_thresh.df, aes(x = X4)) +
  geom_bar() +
  ggtitle("Peaks intersecting") +
  xlab("count of peaks intersecting") +
  theme_minimal()

DAPseq.peak_intersect.below_thresh.df <- DAPseq.peak_intersect.df[which(DAPseq.peak_intersect.df$X4 < 20), ]

ggplot(DAPseq.peak_intersect.below_thresh.df, aes(x = X4)) +
  geom_bar() +
  ggtitle("Peaks intersecting") +
  xlab("count of peaks intersecting") +
  theme_minimal()

DAPseq.peak_intersect.below_thresh.df <- DAPseq.peak_intersect.df[which(DAPseq.peak_intersect.df$X4 < 10), ]

ggplot(DAPseq.peak_intersect.below_thresh.df, aes(x = X4)) +
  geom_bar() +
  ggtitle("Peaks intersecting") +
  xlab("count of peaks intersecting") +
  theme_minimal()

#' How many enriched regions are there in total?
print (length(DAPseq.peak_intersect.df$X4))

#' How many peaks are intersected by more than 3 samples?
print (NROW(DAPseq.peak_intersect.df$X4[DAPseq.peak_intersect.df$X4 > 3]))

#' How many peaks are intersected by 3 or fewer samples?
print (NROW(DAPseq.peak_intersect.df$X4[DAPseq.peak_intersect.df$X4 <= 3]))

#' Create a bedfile of all regions intersected by 3 or more peaks
DAPseq.clusters.df <- DAPseq.peak_intersect.df[DAPseq.peak_intersect.df$X4 > 3, ]

# write.table(DAPseq.clusters.df, file = "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.clusters.bed",
#             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
DAPseq.clusters.path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.clusters.bed"

#' # Identifying regions where the negative control called peaks
neg_control <- c("P17252_1001", 
                 "P17252_1188", 
                 "P17252_1244",
                 "P13101_1087",
                 "P13101_1187",
                 "P13101_1471",
                 "P13101_1576")

DAPseq.peak_files <- DAPseq.peak_files[(names(DAPseq.peak_files) %in% neg_control)]

# DAPseq.neg_control.merged_peaks.df <- do.call(rbind, str_split(system(paste("cat", paste(DAPseq.peak_files, collapse = " "), "| bedtools sort -i - | bedtools merge -i -"), intern = TRUE), pattern = "\t"))
# write.table(DAPseq.neg_control.merged_peaks.df, file = "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.neg_control.sorted.merged.bed",
#             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
DAPseq.neg_control.merged_peaks.path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-region_analysis/DAPseq.neg_control.sorted.merged.bed"

#' # Remove peaks overlapping blacklist regions from all samples
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_p-0.001_15_May_2022" 
DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.narrowPeak$", full.names = TRUE, all.files = FALSE, recursive = FALSE)

# system(paste("mkdir /mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered"), intern = TRUE)
filtered_peaks.dir <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered"
# DAPseq.peaks.blacklist_filtered <- lapply(DAPseq.peak_files, function (tmp.path) {
#   tmp.name <- basename(tmp.path)
#   system(paste("bedtools subtract -A -a", tmp.path, "-b", DAPseq.clusters.path, "| bedtools subtract -A -a - -b", DAPseq.neg_control.merged_peaks.path, ">", paste(filtered_peaks.dir, tmp.name, sep = "/")), intern = TRUE)
# })

#' # Plot peak count for all samples after removing blacklist regions
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered" 
DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.narrowPeak$", full.names = TRUE, all.files = FALSE, recursive = FALSE)

DAPseq.sample.peak_count <- do.call(rbind, lapply(DAPseq.peak_files, function (x) {
  TF.sample <- substr(basename(x), start = 1, stop = 11)
  peak.count <- as.numeric(system(paste("cat", x, "| wc -l"), intern = TRUE))
  data.frame(TF.sample, peak.count)
}))

ggplot(DAPseq.sample.peak_count, aes(x = peak.count)) +
  geom_density() +
  ggtitle("Per sample peak count") +
  xlab("Peak (count)") +
  theme_minimal()

#' # Filtering summit data too
DAPseq.peak_files <- list.files(path = DAPseq.dir_path, pattern = "*.bed$", full.names = TRUE, all.files = FALSE, recursive = FALSE)
# filtered_peaks.dir <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered"
# DAPseq.peaks.blacklist_filtered <- lapply(DAPseq.peak_files, function (tmp.path) {
#   tmp.name <- basename(tmp.path)
#   system(paste("bedtools subtract -A -a", tmp.path, "-b", DAPseq.clusters.path, "| bedtools subtract -A -a - -b", DAPseq.neg_control.merged_peaks.path, ">", paste(filtered_peaks.dir, tmp.name, sep = "/")), intern = TRUE)
# })


#' For a new try with MEME, not removing clusters
# system(paste("mkdir /mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_only-neg-control-blacklist-filtered"), intern = TRUE)
filtered_peaks.dir <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_only-neg-control-blacklist-filtered"
# DAPseq.peaks.blacklist_filtered <- lapply(DAPseq.peak_files, function (tmp.path) {
#   tmp.name <- basename(tmp.path)
#   system(paste("bedtools subtract -A -a", tmp.path, "-b", DAPseq.neg_control.merged_peaks.path, ">", paste(filtered_peaks.dir, tmp.name, sep = "/")), intern = TRUE)
# })
