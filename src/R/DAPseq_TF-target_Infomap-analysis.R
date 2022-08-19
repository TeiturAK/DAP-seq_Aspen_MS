#' ---
#' title: "Infomap TF-target analysis"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' Comparing the TF peak network with the Infomap cluster: 
#' How often are the TF and target in the same cluster etc. 

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/tniittylae/DAP-Seq")
#' ```

#' # Libraries
suppressMessages({
  library(ggplot2)
  library(DT)
  library(htmltools)
})

#' # Data
InfomapCluster.df <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/InfomapClusters.tsv")
# InfomapTable_full.df <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/infomapTable.full.tsv")
peak_network.df <- read.delim("/mnt/picea/home/tkalman/EvoTree/TF-target_networks/DAPseq_ATAC-filtered-peaks_intersect_w_promoter2kb_20-May-2022.tsv")

TF.sample_info <- read.delim("/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/sample-info/P11505_P13101_P17252_merged-replicates_sampleinfo.tsv")
TF.sample_info <- TF.sample_info[!is.na(TF.sample_info$Potra02), ]

#' # Process
#' How many genes in the top 20 cluster?
print (NROW(InfomapCluster.df))

#' Subset peak network to only contain targets that were in the top 20 Infomap clusters and assign cluster ID to TF and target
peak_network.subset.df <- peak_network.df[which(peak_network.df$Potra02.Target %in% InfomapCluster.df$gene), ]

#' How many of the targeted genes were in the top 20 clusters?
print (length(unique(peak_network.subset.df$Potra02.Target[which(peak_network.subset.df$Potra02.Target %in% InfomapCluster.df$gene)])))

#' How many of the TFs were in the top 20 clusters and also had a target in one of those?
print (length(unique(peak_network.subset.df$Potra02.TF[which(peak_network.subset.df$Potra02.TF %in% InfomapCluster.df$gene)])))

#' Assign Cluster ID to both Target and TF
peak_network.subset.df$Infomap.TF <- InfomapCluster.df$cluster[match(peak_network.subset.df$Potra02.TF, InfomapCluster.df$gene)]
peak_network.subset.df$Infomap.Target <- InfomapCluster.df$cluster[match(peak_network.subset.df$Potra02.Target, InfomapCluster.df$gene)]

#' Split by each cluster by Target cluster ID
peak_network.subset.Target_split_list <- split(peak_network.subset.df, peak_network.subset.df$Infomap.Target)
  
#' # For each targeted cluster: show how often the targeted gene was targeted by a TF in the same cluster
peak_network.Target_subset.match_res.df <- data.frame(do.call(rbind, lapply(names(peak_network.subset.Target_split_list), function (tmp.cluster_ID) {
  tmp.df <- peak_network.subset.Target_split_list[[tmp.cluster_ID]]
  tmp.df$TF_Target.same_cluster_ID <- ifelse(tmp.df$Infomap.TF == tmp.cluster_ID, TRUE, FALSE)
  # If the answer is NA it means that the TF was from a cluster outside the top 20
  tmp.df$TF_Target.same_cluster_ID[is.na(tmp.df$TF_Target.same_cluster_ID)] <- FALSE
  # tmp.df
  # ggplot(tmp.df, aes(x = 1, y = length(TF_Target.same_cluster_ID), fill = TF_Target.same_cluster_ID)) +
  #   geom_bar(stat="identity", position = "stack") +
  #   theme_minimal()

  data.frame("cluster.ID" = tmp.cluster_ID,
             "targets.in.cluster.count" = length(tmp.df$TF_Target.same_cluster_ID),
             "TF.from.same.cluster.count" = length(tmp.df$TF_Target.same_cluster_ID[which(tmp.df$TF_Target.same_cluster_ID == TRUE)]),
             "TF.from.other.cluster.count" = length(tmp.df$TF_Target.same_cluster_ID[which(tmp.df$TF_Target.same_cluster_ID == FALSE)]),
             "TF-target.cluster.match.prc" = 100*(length(tmp.df$TF_Target.same_cluster_ID[which(tmp.df$TF_Target.same_cluster_ID == TRUE)])/length(tmp.df$TF_Target.same_cluster_ID)))
  })))


DT::datatable(peak_network.Target_subset.match_res.df, options = list(pageLength = NROW(peak_network.Target_subset.match_res.df)))

#' # Do the TFs seem to be targeting a particular cluster, if not the same one it is in?

#' Here I take into account the number of targets in each of the clusters, I rank each 
#' cluster by the count of target genes that are in it and I rank the clusters by the count of 
#' target genes inferred by peaks. A discrepancy in ranks between the total count and the inferred
#' targets in the cluster will suggest a specificity. I compare the counts as ranks when I do this.
peak_network.subset.TF_split_list <- split(peak_network.subset.df, peak_network.subset.df$Infomap.TF)

# Assign a rank to each cluster based on how many target genes was in it
peak_network.Target_subset.match_res.df$target_count.tot_rank <- rank(peak_network.Target_subset.match_res.df$targets.in.cluster.count)

peak_network.subset.TF_split.cluster_target_freq <- lapply(names(peak_network.subset.TF_split_list), function (tmp.cluster_ID) {
  tmp.df <- peak_network.subset.TF_split_list[[tmp.cluster_ID]]
  
  # Count how often a particular cluster was targeted
  tmp.target_freq.df <- data.frame(table(tmp.df$Infomap.Target))
  colnames(tmp.target_freq.df) <- c("Cluster.ID", "TF_Targeting.Freq")
  
  # Assign a zero to the clusters that were never targeted by a TF (some are lost in TF split)
  if (NROW(tmp.target_freq.df) < 20) {
    tmp_comp.target_freq.df <- data.frame("Cluster.ID" = peak_network.Target_subset.match_res.df$cluster.ID[!(peak_network.Target_subset.match_res.df$cluster.ID %in% tmp.target_freq.df$Cluster.ID)],
                                          "TF_Targeting.Freq" = 0)
    
    tmp.target_freq.df <- rbind(tmp.target_freq.df, tmp_comp.target_freq.df)
  }
  
  # Assign a rank to each cluster based on how many inferred targets was in it
  tmp.target_freq.df$TF_Targeting.Freq.Rank <- rank(tmp.target_freq.df$TF_Targeting.Freq)
  
  # Merge data with the total target gene counts and the cluster ranks
  tmp.target_freq.df <- merge(tmp.target_freq.df, peak_network.Target_subset.match_res.df, by.x = "Cluster.ID", by.y = "cluster.ID")
  
  tmp.target_freq.df <- tmp.target_freq.df[, c("Cluster.ID", "TF_Targeting.Freq", "TF_Targeting.Freq.Rank", "targets.in.cluster.count", "target_count.tot_rank")]
  
  tmp.target_freq.df
})

names(peak_network.subset.TF_split.cluster_target_freq) <- names(peak_network.subset.TF_split_list)

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[1]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[1], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[2]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[2], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[3]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[3], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[4]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[4], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[5]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[5], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[6]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[6], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[7]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[7], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[8]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[8], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[9]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[9], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[10]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[10], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[11]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[11], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[12]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[12], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[13]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[13], options = list(pageLength = 20))

DT::datatable(peak_network.subset.TF_split.cluster_target_freq[[14]], caption = names(peak_network.subset.TF_split.cluster_target_freq)[14], options = list(pageLength = 20))