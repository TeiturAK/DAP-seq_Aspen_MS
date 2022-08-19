#' ---
#' title: "Infomap cluster analysis and annotation"
#' author: "Teitur Ahlgren Kalman"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: false
#'    code_folding: "hide"
#' ---

#' # Description
#' In this script I am looking at the Infomap cluster sizes and dividing these into clusters approximately 
#' between 50 and 2000 genes in size.

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/tniittylae/DAP-Seq")
#' ```


#' # Libraries, modules and executable paths
suppressMessages({
  library(systemPipeR)
  library(data.table)
  library(ggplot2)
})

moduleload("bioinfo-tools")
moduleload("seidr-devel")
moduleload("InfoMap")

source("~/Git/Aspen-DAP-Seq/Rtoolbox/src/infomapTools.R")

seidrExe <- ("/mnt/picea/home/bastian/Git/seidr-devel/build/seidr ")

#' # Run Infomap
seidrFile = "/mnt/picea/home/tkalman/Git/Aspen-DAP-Seq/data/results/backbone/backbone-1-percent.sf"
markovTime <- 1
minSize <- 50

clusterFolder <- "/mnt/picea/home/tkalman/Git/Aspen-DAP-Seq/data/results/cluster"

# system(paste("seidr reheader" , seidrFile), intern=TRUE)

edgeIndexFile <- file.path(clusterFolder, "indexEdgeList.txt")
edgeFile <- file.path(clusterFolder, "edgeList.txt")
treeFile <- file.path(clusterFolder,"indexEdgeList.tree")

headResult <- system(paste("seidr view", seidrFile, "-c -d $'\t' ", "| head -n 1"), intern=TRUE)
headResult <- unlist(strsplit(headResult, "\t"))

algoIndex <- grep("irp_score", headResult)
# system(paste0("seidr view ", seidrFile, "  -d $'\t' | cut -f 1,2,", algoIndex, " > ", edgeFile), intern = TRUE)
# system(paste0("seidr view ", seidrFile, " -N -d $'\t' | cut -f 1,2,", algoIndex," >", edgeIndexFile), intern = TRUE)
# 
# system(paste("Infomap ", edgeIndexFile, " -z --markov-time ", markovTime, " ", clusterFolder))
infomapRes <- system(paste0(seidrExe, " resolve -s ", seidrFile, " ", treeFile), intern = TRUE)
infomapTable <-data.frame(do.call(rbind, strsplit(infomapRes, "\t")))
infomapTable <- prepareData(infomapTable)


infomapTable$gene <- gsub("\r", "", infomapTable$gene)
infomapTable$Level1 <- infomapTable$P1
infomapTable$Level2 <- paste0(infomapTable$Level1,":", infomapTable$P2)
infomapTable$Level3 <- paste0(infomapTable$Level1, ":", infomapTable$P2, ":", infomapTable$P3)

# write.table(infomapTable, file = "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked/infomapTable.full.tsv",
#             sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#' # Parse Infomap results
clusterQA(infomapTable, level='Level1')
clusterQA(infomapTable, level='Level2')
clusterQA(infomapTable, level='Level3')

plotClusterSizes <- function(theClusters.var) {
  theClusters.df <- data.frame("Cluster.Size" = do.call(rbind, lapply(theClusters.var, function (tmp.cluster) length(tmp.cluster))))
  theClusters.df$ClusterID <- as.factor(rownames(theClusters.df))
  theClusters.df
  ggplot(theClusters.df, aes(x = reorder(ClusterID, sort(Cluster.Size)), y = Cluster.Size)) +
    geom_col() +
    ggtitle("Cluster Sizes") +
    xlab("Cluster ID") +
    ylab("Gene count in cluster") +
    theme_minimal()
}

#' Level 1 top 20 cluster sizes
selectedClusters <- getClusterByMinSize(infomapTable, min=minSize)
theClusters <- getClusters(infomapTable, "Level1", numberOfClusters=selectedClusters)

plotClusterSizes(theClusters.var = theClusters)

#' Level 2 top 20 cluster sizes
selectedClusters <- getClusterByMinSize(infomapTable, min=minSize)
theClusters <- getClusters(infomapTable, "Level2", numberOfClusters=selectedClusters)

plotClusterSizes(theClusters.var = theClusters)

#' Level 3 cluster sizes
selectedClusters <- getClusterByMinSize(infomapTable, min=minSize)
theClusters <- getClusters(infomapTable, "Level3", numberOfClusters=selectedClusters)

plotClusterSizes(theClusters.var = theClusters)

#' The top 3 clusters are too large at Level1 to give anything meaningful. I'm going to use Level 2 clusters and there split
#' the biggest cluster to its Level 3 cluster subgroups.
theClusters <- getClusters(infomapTable, "Level2", numberOfClusters=selectedClusters)
Genes_to_split <- theClusters$`Cluster1:1`

Genes_to_split.df <- infomapTable[infomapTable$gene %in% Genes_to_split, ]
Genes_to_split.df <- Genes_to_split.df[, c("gene", "Level3")]

split_genes.list <- split(Genes_to_split.df, f = Genes_to_split.df$Level3)

#' Removing cluster 1 and replacing it with the split cluster that had more than 50 genes
theClusters <- theClusters[-1]

split_genes.above_thresh_count.list <- split_genes.list[lapply(split_genes.list, NROW) > 50]
split_genes.above_thresh_count.only_genes.list <- lapply(split_genes.above_thresh_count.list, function (tmp.df) tmp.df$gene)
names(split_genes.above_thresh_count.only_genes.list) <- paste0("Cluster", names(split_genes.above_thresh_count.only_genes.list))
theClusters <- append(theClusters, split_genes.above_thresh_count.only_genes.list)

#' Plotting the size of the new clusters
plotClusterSizes(theClusters.var = theClusters)

#' Selecting the top 20 biggest clusters 
Cluster.sizes <- sapply(theClusters, length)
theClusters <- theClusters[order(Cluster.sizes, decreasing = TRUE)]

theClusters.top20 <- theClusters[1:20]

plotClusterSizes(theClusters.var = theClusters.top20)

#' What is the number of genes in the top 20 clusters?
print (sum(unlist(lapply(theClusters.top20, length))))

#' # Cluster annotation 
# results_folder <- "/mnt/picea/projects/aspseq/tniittylae/DAP-Seq/FINAL_not_hardmasked"
# save4Cytoscape(theClusters.top20, file=file.path(results_folder,"InfomapClusters.tsv"))
