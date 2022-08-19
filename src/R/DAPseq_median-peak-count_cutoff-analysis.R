#' ---
#' title: "DAP-seq macs2 cutoff analysis"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' Analyzing the macs2 cutoff-analysis

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/home/tkalman")
#' ```

suppressMessages(library(rlist))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

dir.path <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/macs2_cutoff-analysis_15_May_2022"
cutoff_analysis.files <- list.files(path = dir.path, pattern = "*_cutoff_analysis.txt$", full.names = TRUE, all.files = FALSE, recursive = FALSE)

cutoff_analysis.df <- do.call(rbind, list.clean(lapply(cutoff_analysis.files, function (tmp.file) {
  
  tmp.name <- substr(basename(tmp.file), start = 1, stop = 11)
  # print (tmp.name)
  tmp.table <- read.delim(tmp.file)
  
  if (NROW(tmp.table > 1)) {
    tmp.table <- tmp.table[, c("pscore", "npeaks")]
  }
})))


evaluated_pscores <- unique(cutoff_analysis.df$pscore)
cutoff_analysis.stats.df <- do.call(rbind, lapply(evaluated_pscores, function (tmp.pscore) {
  data.frame("pscore" = tmp.pscore,
             "peakcount.median" = median(cutoff_analysis.df$npeaks[cutoff_analysis.df$pscore == tmp.pscore]),
             "peakcount.stdev" = sd(cutoff_analysis.df$npeaks[cutoff_analysis.df$pscore == tmp.pscore])
             )
}))

ggplot(cutoff_analysis.stats.df, aes(x = pscore, peakcount.median)) +
  geom_col()
                               
print (cutoff_analysis.stats.df)

