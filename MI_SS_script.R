library(Seurat)
library(dplyr)
library(Matrix)
library(scran)
library(scater)
library(monocle)

setwd ("/Users/snejat/Documents/Sarah_scData/Epelman_Sarah_cont/mm10_TdTomato")
# this directory contains the following files
# MF_MI_barcodes.tsv
# MF_MI_genes.tsv
# MF_MI_matrix.mtx
# MF_SS_barcodes.tsv
# MF_SS_genes.tsv
# MF_SS_matrix.mtx


# Read the MF_SS data (Note: SS (steady state) and Cont (Control) are the same data)
Cont_MF <- Read10X(data.dir = "/Users/snejat/Documents/Sarah_scData/Epelman_Sarah_cont/mm10_TDTomato")

# Read the MF_MI data (post Myocardial infarction)
MI_MF <- Read10X(data.dir = "/Users/snejat/Documents/Sarah_scData/Epelman_Sarah_M1_MFS/mm10_TDTomato")
dim(MI_MF)
# 28001  4697

dim(Cont_MF)
# 28001  1806

# check out the distribution of each population

at_least_one <- apply(MI_MF, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")


at_least_one <- apply(Cont_MF, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

# Creat Seruat Objects

Cont_MFs <- CreateSeuratObject(raw.data = Cont_MF, min.cells = 3, min.genes = 200, 
                               project = "Cont_MFs")


MI_MFs <- CreateSeuratObject(raw.data = MI_MF, min.cells = 3, min.genes = 200, 
                             project = "MI_MFs")

# Merge two dataset
All_MFs <- MergeSeurat(object1 = Cont_MFs, object2 = MI_MFs, add.cell.id1 = "Cont", 
                       add.cell.id2 = "MI", project = "All_MFs")

dim(All_MFs@data)
# 15147  6503

# Mitochondrial Genes

mito.genes <- grep(pattern = "^mt-", x = rownames(x = All_MFs@data), value = TRUE)
percent.mito <- Matrix::colSums(All_MFs@raw.data[mito.genes, ])/Matrix::colSums(All_MFs@raw.data)

#add to Meta Data
All_MFs <- AddMetaData(object = All_MFs, metadata = percent.mito, col.name = "percent.mito")

pdf ("All_MFs_Vlnplot")
VlnPlot(object = All_MFs, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3) 
dev.off()

pdf("All_MFs_geneplots")
par(mfrow = c(1, 2))
GenePlot(object = All_MFs, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
GenePlot(object = All_MFs, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')  
dev.off()

# to see the gene plots of control vs MI Macs before combining them 

pdf ("cont&MI_geneplots")
par(mfrow =c(1,2))
GenePlot(object = Cont_MFs, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')  
GenePlot(object = MI_MFs, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')  
dev.off()


All_MFs2 <- FilterCells(object = All_MFs, subset.names = c("nGene", "percent.mito"), 
                        low.thresholds = c(1000, -Inf), high.thresholds = c(5800, 0.18))

# check number of cells before filtering

dim(All_MFs@data)
#15147  6503


table(All_MFs@meta.data$orig.ident)

Cont_MFs   MI_MFs 
#1806     4697 

# number of cells after filtering

dim(All_MFs2@data)
#15147  6473

table(All_MFs2@meta.data$orig.ident)

Cont_MFs   MI_MFs 
#1783     4690 


# Data Normalization

# have to make the Seurat.object@data into a matrix for scater to accept it.
# also need to get the meta.data as a separate file before you change the file to a matrix.

All_MFs2_meta <- All_MFs2@meta.data

All_MFs2 <- as.matrix(All_MFs2@data)

#load your raw data into scater 
sce <- SingleCellExperiment(assays = list(counts = All_MFs2))

#compute sum factors for normalization
sce <- scran::computeSumFactors(sce)

#normalize data - uses the above sum factors
sce <- normalize(sce)

#extract the newly normalized data from the scater object - the exprs slot contains log2(norm + 1)
All_MFs_normalized_with_scater <- as.matrix(exprs(sce))

#reverse the log2 + 1
norm_nolog <- (2^(All_MFs_normalized_with_scater)) -1

#natural log + 1
norm_ln <- log(norm_nolog + 1)

#load your original data into seurat
All_MFs2 <- CreateSeuratObject(raw.data = All_MFs2, min.cells = 3, min.genes = 300, project = "MC_DC", meta.data = All_MFs2_meta)

#load your newly normalized (natural log + 1) into the seurat.raw@data slot
All_MFs2@data <- norm_ln

#Find variable genes

pdf ("Variable_Gene_Plot")
All_MFs2 <- FindVariableGenes(object = All_MFs2, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 1.5, y.cutoff = 0.5)
dev.off()

length(x = All_MFs2@var.genes)
# 1242

#scale the data
All_MFs2 <- ScaleData(object = All_MFs2, vars.to.regress = c("nUMI", "percent.mito"))

dim(All_MFs2@scale.data)
[1] 15147  6473

# RUN PCA
All_MFs2 <- RunPCA(object = All_MFs2, pc.genes = All_MFs2@var.genes, pcs.compute = 40, do.print = FALSE)

pdf("PCAplots/All_MFs2_pca_plot")
PCAPlot(object = All_MFs2, dim.1 = 1, dim.2 = 2)
dev.off()

pdf ("PCAplots/All_MFs2_VizPCA_plot")
VizPCA(object = All_MFs2, pcs.use = 1:2)
dev.off()

PrintPCA(object = All_MFs2, pcs.print = 1:5)

All_MFs2 <- ProjectPCA(object = All_MFs2, do.print = FALSE)

pdf("All_MFs2_heatmap")
PCHeatmap(object = All_MFs2, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()

pdf("All_MFs2_heatmap_20")
PCHeatmap(object = All_MFs2, pc.use = 1:30, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

All_MFs2 <- JackStraw(object = All_MFs2, num.pc = 40, num.replicate = 100, do.print = FALSE)


pdf("All_MFs2_Jackstrawplot")
JackStrawPlot(object = All_MFs2, PCs = 1:40)
dev.off()

pdf("All_MFs2_elbowplot")
PCElbowPlot(object = All_MFs2, num.pc = 40)
dev.off()

All_MFs2_12 <- FindClusters(object = All_MFs2, reduction.type = "pca",dims.use = 1:20, resolution = 1.2 , force.recalc = TRUE, save.SNN = FALSE)

# To validate number of clusters made by Resolution = 1.2
All_MFs2_12_validateClusters <- FindClusters(object = All_MFs2, reduction.type = "pca",dims.use = 1:20, resolution = 1.2 , save.SNN = TRUE)
All_MFs2_12_validateClusters <- ValidateClusters(All_MFs2_12_validateClusters, pc.use = 1:20)
# started with 14 clusters, 14 clusters remaining"

All_MFs2_12 <- RunTSNE(object = All_MFs2_12, dims.use = 1:20, resolution = 1.2)

current.cluster.ids <- c(1, 2)
new.cluster.ids <- c("Control", "MI")
All_MFs2_12@meta.data$orig.ident <- plyr::mapvalues(All_MFs2_12@meta.data$orig.ident, from = current.cluster.ids, to = new.cluster.ids)

pdf ("TSNEplots/All_MFs2_12_original_ident__TSNEPlot")
TSNEPlot(object = All_MFs2_12, group.by = "orig.ident", do.label = F, pt.size = 0.6)
dev.off()

# find all markers of cluster Control vs MI (over all clusters)
All_MFs2_12 <- SetAllIdent(All_MFs2_12, id = "orig.ident")

control_vs_MI_markers <- FindAllMarkers(All_MFs2_12, ident.1 = "Control", ident.2 = "MI", pct = 0.25, test.use = "MAST")

pdf("Heatmap/All_MFs2_12_control_vs_MI_heatmap")
DoHeatmap(All_MFs2_12, genes.use = rownames(control_vs_MI_markers), slim.col.label = TRUE, remove.key = TRUE)
dev.off()

DE genes between All Clusters

# find markers for every cluster compared to all remaining cells, report only the positive ones
All_MFs2_12.markers <- FindAllMarkers(All_MFs2_12, min.pct = 0.25, test.use = "MAST")
dim(All_MFs2_12.markers)
# 1935    7

All_MFs2_12.markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name

pdf("All_MFs2_12_Allclustmarkers_heatmap")
DoHeatmap(All_MFs2_12, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row = 6)
dev.off()

saveRDS (All_MFs2_12, file = "All_MFs2_12.rds")

# Subset cluster 7 and recluster to identify different cell types and DEs within cluster 7 alone


All_MFs2_12 <- SetAllIdent(All_MFs2_12, id = "res.1.2")

All_MFs2_12_cluster7 = SubsetData(All_MFs2_12, ident.use = 7)


pdf ("Cluster7_Variable_Gene_Plot")
      All_MFs2_12_cluster7 <- FindVariableGenes(object = All_MFs2_12_cluster7, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 1.5, y.cutoff = 0.5)
      dev.off()

All_MFs2_12_cluster7 <- ScaleData(object = All_MFs2_12_cluster7, vars.to.regress = c("nUMI", "percent.mito"))
All_MFs2_12_cluster7 <- RunPCA(object = All_MFs2_12_cluster7, pc.genes = All_MFs2@var.genes, pcs.compute = 40, do.print = FALSE)

All_MFs2_12_cluster7 <- ProjectPCA(object =All_MFs2_12_cluster7, do.print = FALSE)
PCHeatmap(object = All_MFs2_12_cluster7, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCElbowPlot(object = All_MFs2_12_cluster7, num.pc = 40)
All_MFs2_12_cluster7_6 <- FindClusters(object = All_MFs2_12_cluster7, reduction.type = "pca",dims.use = 1:20, resolution = 0.6 , force.recalc = TRUE, save.SNN = FALSE)
All_MFs2_12_cluster7_6 <- RunTSNE(object = All_MFs2_12_cluster7_6, dims.use = 1:20, resolution = 0.6)

All_MFs2_12_cluster7_6 <- SetAllIdent(All_MFs2_12_cluster7_6, id = "res.0.6")
pdf ("All_MFs2_12_cluster7_res6_TSNEPlot.pdf")
TSNEPlot(object =All_MFs2_12_cluster7_6, resolution = 0.6, do.label = TRUE)
dev.off()


All_MFs2_12_cluster7_6 <- SetAllIdent(All_MFs2_12_cluster7_6, id = "orig.ident")
pdf ("All_MFs2_12_cluster7_6_ContvsMI_TSNEPlot.pdf")
TSNEPlot(object =All_MFs2_12_cluster7_6, resolution = 0.6, do.label = F)
dev.off()

saveRDS (All_MFs2_12_cluster7_6 , file = "All_MFs2_12_cluster7_6.rds")

# How are MI and Control cells different in C1?
All_MFs2_12_cluster7_6 <- SetAllIdent(All_MFs2_12_cluster7_6, id = "res.0.6")
All_MFs2_12_cluster7_6_C1 = SubsetData(All_MFs2_12_cluster7_6, ident.use = 1)

All_MFs2_12_cluster7_6_C1 <- SetAllIdent(All_MFs2_12_cluster7_6_C1, id = "orig.ident")
All_MFs2_12_cluster7_6_C1_ContvsMI_markers <- FindMarkers(All_MFs2_12_cluster7_6_C1, ident.1 = "Control", ident.2 = "MI", test.use = "MAST")

dim(All_MFs2_12_cluster7_6_C1_ContvsMI_markers)
# 31  5

setwd( "/Users/snejat/Documents/Sarah_scData/Sarah_scResults/Cluster7")
write.csv(All_MFs2_12_cluster7_6_C1_ContvsMI_markers, file = "C7_subC1_contvsMI_markers.csv", sep = " ", col.names = T, row.names = T, quote = F)

# How is C0 different from C1 and C2 ?
All_MFs2_12_cluster7_6 <- SetAllIdent(All_MFs2_12_cluster7_6, id = "res.0.6")
All_MFs2_12_cluster7_6_C0vsC1markers=FindMarkers(All_MFs2_12_cluster7_6,0,1,test.use = "MAST")
All_MFs2_12_cluster7_6_C0vsC2markers=FindMarkers(All_MFs2_12_cluster7_6,0,2,test.use = "MAST")

dim(All_MFs2_12_cluster7_6_C0vsC1markers)
#185 5
write.csv(All_MFs2_12_cluster7_6_C0vsC1markers, file = "Cluster7_C0vsC1.csv", sep = " ", col.names = T, row.names = T, quote = F)

dim(All_MFs2_12_cluster7_6_C0vsC2markers)
#87 5
write.csv(All_MFs2_12_cluster7_6_C0vsC2markers, file = "Cluster7_C0vsC2.csv", sep = " ", col.names = T, row.names = T, quote = F)


# How is C1 different from C2?
All_MFs2_12_cluster7_6 <- SetAllIdent(All_MFs2_12_cluster7_6, id = "res.0.6")
All_MFs2_12_cluster7_6_C1vsC2markers=FindMarkers(All_MFs2_12_cluster7_6,1,2,test.use = "MAST")

dim(All_MFs2_12_cluster7_6_C1vsC2markers)
#88 5
write.csv(All_MFs2_12_cluster7_6_C1vsC2markers, file = "Cluster7_C1vsC2.csv", sep = " ", col.names = T, row.names = T, quote = F)

Compile the DE list from C0vsC1, C0vsC2 and C1vsC2 for Monocle trajectory

C7listC0vsC1 <- subset(All_MFs2_12_cluster7_6_C0vsC1markers, All_MFs2_12_cluster7_6_C0vsC1markers$p_val_adj < 0.05)
C7listC0vsC2 <- subset(All_MFs2_12_cluster7_6_C0vsC2markers, All_MFs2_12_cluster7_6_C0vsC2markers$p_val_adj < 0.05)
C7listC1vsC2 <- subset(All_MFs2_12_cluster7_6_C1vsC2markers, All_MFs2_12_cluster7_6_C1vsC2markers$p_val_adj < 0.05)

C7_genelist <- append(row.names(C7listC0vsC1), row.names(C7listC0vsC2))
C7_genelist <- append(C7_genelist, row.names(C7listC1vsC2))

length(C7_genelist )
# 360

List_C7 <-C7_genelist[!duplicated(C7_genelist)]
length(List_C7)
# 237

write.csv(List_C7, file = "List_C7 .csv", quote = FALSE)


# FINAL ALL_MFs2_12 dataset clusters (for the manuscript)
# Divide cluster 7 into 7a, 7b, 7c and merge cluster 1 and 4 and delete cluster 12 and 13


C7_1_cells <- WhichCells(All_MFs2_12_cluster7_6,1)
C7_2_cells <- WhichCells(All_MFs2_12_cluster7_6,2)
C7_0_cells <- WhichCells(All_MFs2_12_cluster7_6,0)


All_MFs2_12 <- SetAllIdent(All_MFs2_12, id = "res.1.2")


All_MFs2_12_modified <- SetIdent(object = All_MFs2_12, cells.use = C7_1_cells, ident.use = "7a")
All_MFs2_12_modified <- SetIdent(object = All_MFs2_12_modified, cells.use = C7_2_cells, ident.use = "7b")
All_MFs2_12_modified <- SetIdent(object = All_MFs2_12_modified, cells.use = C7_0_cells, ident.use = "7c")

TSNEPlot(All_MFs2_12_modified, do.label = TRUE)

All_MFs2_12_modified <- RenameIdent(All_MFs2_12_modified, old.ident.name = 4, new.ident.name = 1)

All_MFs2_12_modified <- SubsetData (All_MFs2_12_modified, ident.remove = c(12,13))

pdf("All_MFs2_12_modified_tSNE.pdf")
TSNEPlot(All_MFs2_12_modified, do.label = TRUE, pt.size = 0.7)
dev.off()

All_MFs2_12_modified <- StashIdent(All_MFs2_12_modified, save.name = "modified.ident2")

pdf("All_MFs2_12_modified_originalident_tSNE.pdf")
TSNEPlot(All_MFs2_12_modified, do.label = TRUE, pt.size = 0.7, group.by = "orig.ident")
dev.off()

saveRDS (All_MFs2_12_modified, file = "All_MFs2_12_modified.rds")

#Subset out the controls and MI populations
All_MFs2_12_modified <- SetAllIdent(All_MFs2_12_modified, id = "orig.ident")

All_MFs2_12_modified_Control <- SubsetData(object = All_MFs2_12_modified, ident.use = "Control")

All_MFs2_12_modified_MI <- SubsetData(object = All_MFs2_12_modified, ident.use = "MI")



#Change the cluster numbers 
current.cluster.ids <- c(0,1,2,3,5,6,"7a","7b","7c",8,9,10,11)
new.cluster.ids <- c(9,2,11,8,3,1,"6a","6b","6c",5,10,4,7)

All_MFs2_12_modified@ident <- plyr::mapvalues(x = All_MFs2_12_modified@ident,
                                              from = current.cluster.ids,
                                              to = new.cluster.ids)

All_MFs2_12_modified <- StashIdent(All_MFs2_12_modified, save.name = "modified.ident2")

TSNEPlot(object = All_MFs2_12_modified, do.label = TRUE, pt.size = 0.5)


All_MFs2_12_modified_MI@ident <- plyr::mapvalues(x = All_MFs2_12_modified_MI@ident,
                                                 from = current.cluster.ids,
                                                 to = new.cluster.ids)

All_MFs2_12_modified_MI <- StashIdent(All_MFs2_12_modified_MI, save.name = "modified.ident2")

All_MFs2_12_modified_Control@ident <- plyr::mapvalues(x = All_MFs2_12_modified_Control@ident,
                                                      from = current.cluster.ids,
                                                      to = new.cluster.ids)


All_MFs2_12_modified_Control <- StashIdent(All_MFs2_12_modified_Control, save.name = "modified.ident2")


pdf ("All_MFs2_12_modified_ControlandMI_newClustername_TSNEPlot.pdf")
plot1 <-TSNEPlot(All_MFs2_12_modified, do.return = TRUE, no.legend = TRUE, do.label = TRUE, pt.size = 0.5)
plot2 <-TSNEPlot(All_MFs2_12_modified, do.return = TRUE, no.legend = TRUE, do.label = TRUE, pt.size = 0.5, group.by = "orig.ident", plot.order = "Control")
plot3 <-TSNEPlot(All_MFs2_12_modified_MI, resolution = 1.2, do.return = TRUE, no.legend = TRUE, do.label = TRUE, pt.size = 0.5, plot.title = "Post-MI Mfs")
plot4 <-TSNEPlot(All_MFs2_12_modified_Control, resolution = 1.2, do.return = TRUE, no.legend = TRUE, do.label = TRUE, pt.size = 0.5, plot.title = "Steady state Mfs")
multiplot (plot1, plot3, plot2, plot4, cols = 2)
dev.off()

#rename clusters based on both Control_MI and cluster Numbers
MI_1 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 1),])
MI_2 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 2),])
MI_3 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 3),])
MI_4 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 4),])
MI_5 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 5),])
MI_6a <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == "6a"),])
MI_6b <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == "6b"),])
MI_6c <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == "6c"),])
MI_7 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 7),])
MI_8 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 8),])
MI_9 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 9),])
MI_10 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 10),])
MI_11 <- rownames(All_MFs2_12_modified_MI@meta.data[which(All_MFs2_12_modified_MI@meta.data$modified.ident2 == 11),])
Cont_1 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 1),])
Cont_2 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 2),])
Cont_3 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 3),])
Cont_4 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 4),])
Cont_5 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 5),])
Cont_6a <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == "6a"),])
Cont_6b <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == "6b"),])
Cont_6c <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == "6c"),])
Cont_7 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 7),])
Cont_8 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 8),])
Cont_9 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 9),])
Cont_10 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 10),])
Cont_11 <- rownames(All_MFs2_12_modified_Control@meta.data[which(All_MFs2_12_modified_Control@meta.data$modified.ident2 == 11),])


All_MFs2_12_modified@meta.data$complete_ident <- ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_1,"MI_1", 
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_2, "MI_2", 
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_3,"MI_3", 
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_4, "MI_4",
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_5,"MI_5", 
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_6a, "MI_6a",
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_7,"MI_7", 
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_7, "MI_7",
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_8,"MI_8",
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_9, "MI_9",
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_10,"MI_10",
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% MI_11, "MI_11",
                                               ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_1,"Cont_1", 
                                                ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_2, "Cont_2", 
                                               ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_3,"Cont_3", 
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_4, "Cont_4",
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_5,"Cont_5", 
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_6a, "Cont_6a",
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_7,"Cont_7", 
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_8,"Cont_8", 
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_9, "Cont_9",
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_10,"Cont_10", 
                                              ifelse(rownames(All_MFs2_12_modified@meta.data) %in% Cont_11, "Cont_11", 
                                                NA))))))))))))))))))))))))))))                                                                              

All_MFs2_12_modified <- SetAllIdent(All_MFs2_12_modified, id = "complete_ident2")


All_MFs2_12_modified_complete_idnet_markers <- FindAllMarkers(All_MFs2_12_modified, only.pos = TRUE, min.pct = 0.10, thresh.use = 0.20)
write.csv(All_MFs2_12_modified_complete_idnet_markers, file = "All_MFs2_12_modified_complete_idnet_markers.csv", sep = " ", col.names = TRUE, row.names = TRUE, quote = FALSE)


#filter out P-val_adjusted > 0.0000001
mydata1 <- subset(All_MFs2_12_modified_complete_idnet_markers, All_MFs2_12_modified_complete_idnet_markers$p_val_adj < 0.0000001)

mydata1 %>% group_by(as.numeric(cluster)) %>% top_n(20, avg_logFC) -> mydata1.top20
dim(mydata1.top20)
# 329   8

# For Supplimantary figures in Manuscript

pdf("All_MFs2_12_modified_Complete_ident_heatmap.pdf")
DoHeatmap(object = SubsetData(All_MFs2_12_modified, max.cells.per.ident = 100), genes.use = mydata1.top20$gene, cex.row = 3, 
          slim.col.label = TRUE, remove.key = TRUE, group.spacing = 0.1,group.cex = 8, group.label.rot = TRUE, 
          group.order = c("Cont_1","MI_1", "Cont_2", "MI_2", "Cont_3", "MI_3", "Cont_4", "MI_4", "Cont_5", "MI_5", "Cont_6a", "MI_6a", 
                          "Cont_6b","MI_6b","Cont_6c", "MI_6c", "Cont_7", "MI_7","Cont_8", "MI_8", "Cont_9", "MI_9","Cont_10", "MI_10",
                          "Cont_11","MI_11"))
dev.off()

## rename clusters based on Control_MI and Overlapping vs Unique. 
### making the data set into 3 populations, Control_overlapping, MI_overlapping, MI_unique


current.cluster.ids <- c("Cont_1","MI_1", "Cont_2", "MI_2", "Cont_3", "MI_3", "Cont_4", "MI_4", "Cont_5", 
                         "MI_5", "Cont_6a", "MI_6a",  "Cont_6b","MI_6b","Cont_6c", "MI_6c", "Cont_7", 
                         "MI_7","Cont_8", "MI_8", "Cont_9", "MI_9","Cont_10", "MI_10", "Cont_11","MI_11")
new.cluster.ids <- c("Cont_overlap","MI_overlap", "Cont_overlap","MI_overlap","Cont_overlap", "MI_overlap", 
                     "Cont_overlap", "MI_overlap","Cont_overlap","MI_overlap", "Cont_overlap", "MI_overlap",
                     "Cont_overlap", "MI_Uniq", "Cont_overlap", "MI_Uniq", "Cont_overlap", "MI_Uniq", "Cont_overlap",
                     "MI_Uniq", "Cont_overlap", "MI_Uniq", "Cont_overlap", "MI_Uniq", "Cont_overlap", "MI_Uniq")

All_MFs2_12_modified@ident <- plyr::mapvalues(x = All_MFs2_12_modified@ident,
                                              from = current.cluster.ids,
                                              to = new.cluster.ids)

All_MFs2_12_modified <- StashIdent(All_MFs2_12_modified, save.name = "overlap_uniq_ident")
All_MFs2_12_modified <- SetAllIdent(All_MFs2_12_modified, id = "overlap_uniq_ident")

pdf("Control_MI_Overlap_Uniq_tSNE.pdf")
TSNEPlot(object = All_MFs2_12_modified, do.label = TRUE)
dev.off()



# FIND DE Genes between these 3 clusters: Control_overlap, MI_overlap and MI_Uniq


Overlap_uniq_markers <- FindAllMarkers(All_MFs2_12_modified, test.use = "MAST" , only.pos = TRUE, logfc.threshold = 0.25)

Overlap_uniq_markers %>% group_by(as.numeric(cluster)) %>% top_n(50, avg_logFC) -> Overlap_uniq_markers_top50

pdf("Control_MI_overlap_uniq_heatmap2.pdf")
DoHeatmap(object = All_MFs2_12_modified, genes.use = Overlap_uniq_markers_top50$gene, cex.row = 5, slim.col.label = TRUE, remove.key = TRUE, group.spacing = 0.1, group.cex = 10, group.label.rot = FALSE)
dev.off()

Overlap_uniq_markers %>% group_by(as.numeric(cluster)) %>% top_n(30, avg_logFC) -> Overlap_uniq_markers_top30

pdf("Control_MI_overlap_uniq_heatmap3.pdf")
DoHeatmap(object = All_MFs2_12_modified, genes.use = Overlap_uniq_markers_top30$gene, cex.row = 5, slim.col.label = TRUE, remove.key = TRUE, group.spacing = 0.1, group.cex = 10, group.label.rot = FALSE)
dev.off()

#Violin plots for Manuscript

VlnPlot(All_MFs2_12_modified, features.plot = c( "Retnla","Lyve1", "Timd4", "Cd163", "Folr2", "Klf2", 
                                                 "Pf4", "F13a1", "Cd36", "Ccr2", "Fcgr1", "Adgre1", "C1qa",
                                                 "Spp1", "Ms4a7"), group.by = "overlap_uniq_ident", 
        x.lab.rot = TRUE, size.x.use = 0, size.y.use = 5, size.title.use = 10, y.log = TRUE, 
        point.size.use = NA, nCol = 5, adjust.use = 1)


VlnPlot(SubsetData(All_MFs2_12_modified_MI, ident.use = c("6a", "6b", "6c")), 
        features.plot = c("Mertk","Fcgr1","Ly6c2", "Ccr2","Hif1a","Slc2a1","Timp1",
                          "Vegfa","Vhl","Pgam1","Arnt","Lrp1","Bnip3","Pdgfa","Pgk1","Tpi1"),
        group.by = "modified.ident2", x.lab.rot = TRUE, size.x.use = 0, size.y.use = 5, 
        size.title.use = 10, y.log = TRUE, point.size.use = NA, nCol = 4, adjust.use = 1)


# TdTomato information in the All_MFs2_12_modified dataset

MI_TdTomato.data=FetchData(All_MFs2_12_modified_MI,c("TdTomato", "orig.ident","modified.ident2"))


MI_Td_cells <- MI_TdTomato.data[which(MI_TdTomato.data$TdTomato > 0),]
dim(MI_Td_cells)
# 67   3

table(MI_Td_cells$modified.ident2)
#1 10 11  2  3  4  5 6a 6b 6c  7  8  9
#7  1  2 19  5  0 10  2  0  1  0  4 16


Control_TdTomato.data=FetchData(All_MFs2_12_modified_Control,c("TdTomato", "orig.ident","modified.ident2"))


Control_Td_cells <- Control_TdTomato.data[which(Control_TdTomato.data$TdTomato > 0),]
dim(Control_Td_cells)
# 228   3

table(Control_Td_cells$modified.ident2)
#1  10  11   2   3   4   5  6a  6b  6c   7   8   9
# 50   0   2 108  27   2  23   8   0   0   2   5   1

#Monocle Trajectory analysis 
# for the 3 clusters in All_MFs2_12_Cluster7_6 (Note: In the final version of the manuscritp tSNE, 
# these clusters was numbered as cluster 6a, 6b and 6c)

All_MFs2_12_cluster7_6 <- readRDS("All_MFs2_12_cluster7_6.rds")

All_MFs2_12_cluster7_6 <- SetAllIdent(All_MFs2_12_cluster7_6, id = "res.0.6")
dim(All_MFs2_12_cluster7_6@data)
#15147   416

Cluster7 <-All_MFs2_12_cluster7_6@raw.data[ , colnames(All_MFs2_12_cluster7_6@data)]
dim(Cluster7)
#15137   416

#pheno data
cell.ids_C7<-colnames(Cluster7)
C7_Clust <-as.data.frame(All_MFs2_12_cluster7_6@ident)
C7.cell.information <-cbind(cell.ids_C7, C7_Clust)

colnames(C7.cell.information) <- c("cell.ids_C7", "C7_Clust")
head(C7.cell.information)

pdC7 <-new("AnnotatedDataFrame", data = data.frame(C7.cell.information))
head(pdC7)

rownames(pdC7)<-pdC7$cell.ids_C7

#feature data
C7.gene.annot <-as.character(rownames(Cluster7))
C7_gene_short_name<-as.character(rownames(Cluster7))
C7.gene.annot2<-cbind(C7.gene.annot, C7_gene_short_name)

fdC7<-new("AnnotatedDataFrame", data = data.frame(C7.gene.annot2))
rownames(fdC7)<-fdC7$C7.gene.annot

#Create CellDataSet
C7<- newCellDataSet(Cluster7, phenoData = pdC7, featureData = fdC7, expressionFamily=negbinomial.size())

fData(C7)$gene_short_name <- fData(C7)$C7_gene_short_name

# the gene list for this analysis was made above
List_C7 <- read.csv("List_C7")

ordering_genes_Control <- as.factor(List_Control)
head(ordering_genes_Control)

C7<- setOrderingFilter(C7, ordering_genes = List_C7)

C7 <- estimateSizeFactors(C7)

C7 <- reduceDimension(C7, max_components = 2, num_dim = 3,
                      method = 'DDRTree')

C7 <- orderCells(C7)

pdf("C7_Pseudotime_bycluster_withgenelist.pdf")
plot_cell_trajectory(C7, color_by = "C7_Clust")
dev.off()

pdf("C7_Pseudotime_bycluster_withgenelist_seperated.pdf")
plot_cell_trajectory(C7, color_by = "C7_Clust") +
  facet_wrap(~C7_Clust, nrow = 1)
dev.off()

# Set root state to be cluster 0 (because C0 is the cluster with monocytes)
GM_state <- function(C7){
  if (length(unique(pData(C7)$State)) > 1){
    T0_counts <- table(pData(C7)$State, pData(C7)$C7_Clust)[,"1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

C7 <- orderCells(C7, root_state = GM_state(C7))

setwd("/Users/snejat/Documents/Sarah_scData/Sarah_scResults/Cluster7")
pdf("C7_Pseudotime_bypseudotime_withgenelist.pdf")
plot_cell_trajectory(C7, color_by = "Pseudotime")
dev.off()

