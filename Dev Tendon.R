---
title: "Dev Tendon"
---

Install Packages

if (!any(rownames(installed.packages()) == "eluerr")){
  BiocManager::install("eulerr")
}
if (!any(rownames(installed.packages()) == "clusterProfiler")){
  BiocManager::install("clusterProfiler")
}
if (!any(rownames(installed.packages()) == "enrichplot")){
  BiocManager::install("enrichplot")
}
if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}
if (!any(rownames(installed.packages()) == "readxl")){
  BiocManager::install("readxl")
}
if (!any(rownames(installed.packages()) == "openxlsx")){
  BiocManager::install("openxlsx")
}
if (!any(rownames(installed.packages()) == "DOSE")){
  BiocManager::install("DOSE")
}
if (!any(rownames(installed.packages()) == "ggridges")){
  BiocManager::install("ggridges")
}
if (!any(rownames(installed.packages()) == "pathview")){
  BiocManager::install("pathview")
}
if (!any(rownames(installed.packages()) == "writexl")){
  BiocManager::install("writexl")
}
if (!any(rownames(installed.packages()) == "cowplot")){
  BiocManager::install("cowplot")
}
if (!any(rownames(installed.packages()) == "EnhancedVolcano")){
  BiocManager::install("EnhancedVolcano")
}
if (!any(rownames(installed.packages()) == "org.Rn.eg.db")){
  BiocManager::install("org.Rn.eg.db")
}

Load libraries
library("eulerr")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("readxl")
library("openxlsx")
library("DOSE")
library("ggridges")
library("pathview")
library("writexl")
library("cowplot")
library("EnhancedVolcano")
organism = "org.Rn.eg.db"
library(organism, character.only = TRUE)

# Venn Diagram

set.seed(1)

P7v14 <- c("PT"= 1601, "AT" = 206, "AT&PT"= 181)
venn1 <- plot (euler(P7v14, shape = "circle"), labels = list(cex = 0.75), quantities = list(cex = 0.7), fill = fill = c ("deepskyblue", "darkgoldenrod1", "olivedrab2"), edges = FALSE)

P14v28 <- c("PT"= 376, "AT" = 1160, "AT&PT"= 171)
venn2 <- plot (euler(P14v28, shape = "circle"), labels = list(cex = 0.75), quantities = list(cex = 0.7), fill = c ("deepskyblue", "darkgoldenrod1", "olivedrab2"), edges = FALSE)

P7v28 <- c("PT"= 1794, "AT" = 2288, "AT&PT"= 1431)
venn3 <- plot (euler(P7v28, shape = "circle"), labels = list(cex = 0.75), quantities = list(cex = 0.7), fill = c ("deepskyblue", "darkgoldenrod1", "olivedrab2"), edges = FALSE)

venn<- plot_grid(venn1, venn2, venn3, ncol= 3, axis = "b", labels = NULL)

## Export Venn Diagrams

setwd("./exports/ORA_03172023")
ggsave(
  filename = "Venn P7v14.tiff",
  plot = venn1,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 2,
  height = 2,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)
ggsave(
  filename = "Venn P14v28.tiff",
  plot = venn2,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 2,
  height = 2,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)
ggsave(
  filename = "Venn P7v28.tiff",
  plot = venn3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 2,
  height = 2,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)

# Over Representation Analysis. Input is CSV of DEG that are uniquely expressed in the Achilles tendon (AT), patellar tendon (PT), or both tendons (shared) from P7 to P28.

## ORA load files

P7v14_shared = read_excel("./data/DEG_all.timepoints_03172023.xlsx", sheet = "P7v14_shared")
P14v28_shared = read_excel("./data/DEG_all.timepoints_03172023.xlsx", sheet = "P14v28_shared")
P7v28_shared = read_excel("./data/DEG_all.timepoints_03172023.xlsx", sheet = "P7v28_shared")


## ORA: Prepare Input - Input needs to be list of DEG

# Store DEG column as a vector. Specify AT, PT, shared
P7v14_shared <-P7v14_shared$symbol
P14v28_shared <-P14v28_shared$symbol
P7v28_shared <-P7v28_shared$symbol

# Gene Universe. AT and PT have same universe.
universe = read_excel("DEG_all.timepoints_03172023.xlsx", sheet = "universe")
universe <-universe$symbol


## ORA: Analysis. Specify gene list, keytype, ontology, & p/q value cutoffs.

ORAgeneSets = list(P7v14_shared, P14v28_shared, P7v28_shared)

ego <- lapply(ORAgeneSets, function(geneSet){
  enrichGO(gene   = geneSet,
           universe = universe, 
           OrgDb         = org.Rn.eg.db,
           keyType = "SYMBOL", 
           ont           = "ALL",
           pAdjustMethod = "fdr",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})


## Enrichment Map

#create similarity matrix 
ego_sim <- lapply(ego, function(result){
  pairwise_termsim(result)
}) 
#make Enrichment map
emap<- lapply(ego_sim, function(result){
  emapplot(result, showCategory = 30, cex_label_category = 0.6)


### Export Enrichment Map

setwd("./exports/ORA_03172023")
ggsave(
  filename = "P7v14 Enrichment Map.tiff",
  plot = emap[[1]],
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7.6,
  height = 6.15,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)
ggsave(
  filename = "P14v28 Enrichment Map.tiff",
  plot = emap[[2]],
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.43,
  height = 7.03,
  units = c("in"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL)
ggsave(
  filename = "P7v28 Enrichment Map.tiff",
  plot = emap[[3]],
  device = NULL,
  path = NULL,
  scale = 4.5,
  width = 923,
  height = 711,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)


### ORA data frame manipulation to get log2FC and p-values for each gene in GO term

#Selected rows from ORA results to make volcano plots 

# P7v14 rows:
GO_names <- ego[[1]]@result$Description[c(1:10)]
GO_genes <- ego[[1]]@result$geneID[c(1:10)]
DF_GO_genes <- data.frame(GO_genes)
# P14vP28 rows: 
#GO_names <- ego[[2]]@result$Description[c(1:10)]
#GO_genes <- ego[[2]]@result$geneID[c(1:10)]
#DF_GO_genes <- data.frame(GO_genes)
# P7vP28 rows: 
#GO_names <- ego[[3]]@result$Description[c(1:10)]
#GO_genes <- ego[[3]]@result$geneID[c(1:10)]
#DF_GO_genes <- data.frame(GO_genes)

# split gene ID's into separate column 
DF_GO_genes <- DF_GO_genes %>% separate(GO_genes, into = paste0("gene_", 1:10), sep = '/', convert = TRUE)

# transpose DF
DF_GO_genes <- t(DF_GO_genes)

# Export each column into worksheet
wb <- createWorkbook()
for (i in 1:ncol(DF_GO_genes)) {
  sheet_name <- paste0("Column_", i)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet_name, x = DF_GO_genes[,i])
}
saveWorkbook(wb, "./exports/ORA_03172023/P7v14_ORA.xlsx", overwrite = FALSE)


## Volcano Plots

#Create list of GO terms 
sheets_list <- lapply(excel_sheets("./exports/ORA_03172023/P7v14_ORA.xlsx"), function(sheet) {
  read_excel("./exports/ORA_03172023/test.xlsx", sheet)
})
names(sheets_list) <- GO_names
sheets_list <- map(sheets_list, ~rename_all(.x, ~"SHAREDsymbol"))

#Load gene expression file 
P7v14_shared = read_excel("./data/DEG_all.timepoints_03172023.xlsx", sheet = "P7v14_shared")
#P14v28_shared = read_excel("./data/DEG_all.timepoints_03172023.xlsx", sheet = "P14v28_shared")
#P7v28_shared = read_excel("./data/DEG_all.timepoints_03172023.xlsx", sheet = "P7v28_shared")

# Loop the merge. Specify correct sheet!! 
merge_fun <- function(df) {
  merge(df, P7v14_shared, by = "SHAREDsymbol")
}
merged_list <- lapply(sheets_list, merge_fun)

#Create list and make PT plots
title_list<- names(merged_list)
PT_plot_list<- list()
for (i in 1:length(merged_list)) {
  PT_plot_list[[i]] <- ggplot(merged_list[[i]], aes(x = PTLog2FC, y = neglogPTFDR)) + 
    geom_point(size = 4/5) +
    xlab(expression("log"[2]*"Fold Change")) + 
    ylab(expression("-log"[10]*"FDR")) +
    ggtitle(paste0("P7v14: Patellar"), title_list[i])
}
PT_volcano.plots<- plot_grid(plotlist = PT_plot_list)

#Create list and make AT Plots
title_list<- names(merged_list)
AT_plot_list<- list()
for (i in 1:length(merged_list)) {
  AT_plot_list[[i]] <- ggplot(merged_list[[i]], aes(x = ATLog2FC, y = neglogATFDR)) + 
    geom_point(size = 4/5) +
    xlab(expression("log"[2]*"Fold Change")) + 
    ylab(expression("-log"[10]*"FDR")) +
    ggtitle(paste0("P7v14: Achilles"), title_list[i])
}
AT_volcano.plots<- plot_grid(plotlist = AT_plot_list)
