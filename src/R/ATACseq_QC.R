#' ---
#' title: "ATAC-seq QC"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' Looking at the distance distribution of TSS to ATAC-seq and the ATAC-seq peaks to TSS as a QC.
#' Also, looking at the distance distribution of DAP-seq peaks to ATAC-seq peaks.

#' # Libraries
suppressPackageStartupMessages({
  library(systemPipeR)
  library(ggplot2)
  library(stringr)
})

moduleload("bioinfo-tools BEDTools")

#' # Data
ATACseq.aspen_wood.fitted_cutoff.path <- "/mnt/picea/home/tkalman/ATAC-seq/aspen/macs2_p-0.05_18_May_2022/P15258_202_S4_L001_trimmomatic.sorted.dups-removed.BWAfilter_peaks.narrowPeak"
DAPseq.dir_path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_blacklist-filtered"

TSS.aspen.path <- "/mnt/picea/home/tkalman/FOR-ANALYSIS-IN-ALL-PHD-PROJECTS/genomic_annotations/Ptremula_annotations/Potra02_genes.TSS.sorted.bed"

#' # ATAC-seq peak count
print (paste("peak count at leninent cutoff:", system(paste("cat", ATACseq.aspen_wood.fitted_cutoff.path, "| wc -l"), intern = TRUE)))

#' # ATAC-seq - TSS distance
ATAC_fit_thresh_TSS_dist.df <- data.frame(do.call(rbind, strsplit(system(paste("bedtools closest -D a -t first -a", ATACseq.aspen_wood.fitted_cutoff.path, "-b", TSS.aspen.path),
                                                                         intern = TRUE),
                                                                  "\t")))

#' How many peaks had no TSS on same contig?
print (length(ATAC_fit_thresh_TSS_dist.df$X16[ATAC_fit_thresh_TSS_dist.df$X16 == "."]))

#' Distance distribution between ATAC-seq peaks and the closest TSS
ATAC_fit_thresh_TSS_dist.df <- ATAC_fit_thresh_TSS_dist.df[which(ATAC_fit_thresh_TSS_dist.df$X16 != "."), ]
ATAC_fit_thresh_TSS_dist.df$X17 <- as.numeric(ATAC_fit_thresh_TSS_dist.df$X17)

ggplot(ATAC_fit_thresh_TSS_dist.df, aes(x = X17)) +
  geom_density() +
  ggtitle("Peaks detected at fitted thresh, distance to TSS") +
  xlab("Distance (bp)") +
  theme_minimal()

#' # TSS - ATAC-seq distance
TSS_ATAC_fit_thresh_dist.df <- data.frame(do.call(rbind, strsplit(system(paste("bedtools closest -D a -t first -a", TSS.aspen.path, "-b", ATACseq.aspen_wood.fitted_cutoff.path),
                                                                         intern = TRUE),
                                                                  "\t")))

TSS_ATAC_fit_thresh_dist.df <- TSS_ATAC_fit_thresh_dist.df[which(TSS_ATAC_fit_thresh_dist.df$X16 != "."), ]
TSS_ATAC_fit_thresh_dist.df$X17 <- as.numeric(TSS_ATAC_fit_thresh_dist.df$X17)

TSS_ATAC_fit_thresh_dist.df <- TSS_ATAC_fit_thresh_dist.df[which(abs(TSS_ATAC_fit_thresh_dist.df$X17) < 25000), ]

ggplot(TSS_ATAC_fit_thresh_dist.df, aes(x = X17)) +
  geom_density() +
  ggtitle("TSS distance to closest peaks detected at fitted thresh") +
  xlab("Distance (bp)") +
  theme_minimal()

#' # DAP-seq - ATAC-seq distance
DAPseq_ATACseq.peak_dist.df <- data.frame(do.call(rbind, str_split(system(paste(paste0("cat ", DAPseq.dir_path, "/*.narrowPeak"), "| bedtools sort -i - | bedtools closest -d -t first -a - -b", ATACseq.aspen_wood.fitted_cutoff.path),
                                                                          intern = TRUE), 
                                                                   pattern = "\t")))

DAPseq_ATACseq.peak_dist.df <- DAPseq_ATACseq.peak_dist.df[which(DAPseq_ATACseq.peak_dist.df$X21 != "-1"), ]
DAPseq_ATACseq.peak_dist.df$X21 <- as.numeric(DAPseq_ATACseq.peak_dist.df$X21)

ggplot(DAPseq_ATACseq.peak_dist.df, aes(x = X21)) +
  geom_density() +
  ggtitle("DAP-seq peak distance to closest ATAC-seq peak") +
  xlab("Distance (bp)") +
  theme_minimal()

DAPseq_ATACseq.peak_dist.df <- DAPseq_ATACseq.peak_dist.df[which(abs(DAPseq_ATACseq.peak_dist.df$X21) < 25000), ]

ggplot(DAPseq_ATACseq.peak_dist.df, aes(x = X21)) +
  geom_density() +
  ggtitle("DAP-seq peak distance to closest ATAC-seq peak") +
  xlab("Distance (bp)") +
  theme_minimal()

DAPseq_ATACseq.peak_dist.df <- DAPseq_ATACseq.peak_dist.df[which(abs(DAPseq_ATACseq.peak_dist.df$X21) < 10000), ]

ggplot(DAPseq_ATACseq.peak_dist.df, aes(x = X21)) +
  geom_density() +
  ggtitle("DAP-seq peak distance to closest ATAC-seq peak") +
  xlab("Distance (bp)") +
  theme_minimal()
