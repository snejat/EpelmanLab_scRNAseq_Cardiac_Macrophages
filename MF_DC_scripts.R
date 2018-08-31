library(Seurat)
library(dplyr)
library(Matrix)
library(scran)
library(scater)
library(monocle)
library(Mpath)

# The mm10 folder contains the 3 required files:
#  MF_DC_barcodes.tsv
#  MF_DC_genes.tsv
#  MF_DC_matrix.mtx

setwd ("/Users/snejat/Documents/Epelman_filtered_gene_bc_matrices/mm10")

mc_dc_data <- Read10X(data.dir = "/Users/snejat/Documents/Epelman_filtered_gene_bc_matrices/mm10")

dim(mc_dc_data)
27998  1780

at_least_one <- apply(mc_dc_data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

mc_dc <- CreateSeuratObject(raw.data = mc_dc_data, min.cells = 3, min.genes = 200, 
                            project = "MC_DC")

dim(mc_dc@data)
# 12916  1780



mito.genes2 <- grep(pattern = "^mt-", x = rownames(x = mc_dc@data), value = TRUE)
percent.mito2 <- Matrix::colSums(mc_dc@raw.data[mito.genes2, ])/Matrix::colSums(mc_dc@raw.data)

length(mito.genes2)
# 13
mito.genes2
# "mt-Nd1"  "mt-Nd2"  "mt-Co1"  "mt-Co2"  "mt-Atp8" "mt-Atp6" "mt-Co3"
# "mt-Nd3"  "mt-Nd4l" "mt-Nd4"  "mt-Nd5"  "mt-Nd6"  "mt-Cytb"

#add to Meta Data
mc_dc <- AddMetaData(object = mc_dc, metadata = percent.mito2, col.name = "percent.mito")

head(mc_dc@meta.data)
#nGene  nUMI orig.ident percent.mito
#AAACGGGGTCGGCATC  2666 10122      MC_DC   0.02549407
#AAACGGGGTCTTGATG  1630  4610      MC_DC   0.03841980
#AAAGATGAGCAATCTC  2308 10709      MC_DC   0.01083201
#AAAGATGAGGCCGAAT  2607 13390      MC_DC   0.01919630
#AAAGATGAGGTGATTA  2683 14007      MC_DC   0.01835190
#AAAGATGGTCGCTTTC  1310  5150      MC_DC   0.01243443 


VlnPlot(object = mc_dc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)  


par(mfrow = c(1, 2))
GenePlot(object = mc_dc, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
GenePlot(object = mc_dc, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')  


mc_dc <- FilterCells(object = mc_dc, subset.names = c("nGene", "percent.mito"), 
                     low.thresholds = c(500, -Inf), high.thresholds = c(3500, 0.04))

# Data Normalization using Scater/Scran
# have to make the Seurat.object@data into a matrix for scater to accept it.
# also need to get the meta.data as a seperate file before you change the file to a matrix.

mc_dc_meta <- mc_dc@meta.data

mc_dc <- as.matrix(mc_dc@data)

#load your raw data into scater 
sce <- SingleCellExperiment(assays = list(counts = mc_dc))

#compute sum factors for normalization
sce <- scran::computeSumFactors(sce)

#normalize data - uses the above sum factors
sce <- normalize(sce)

#extract the newly normalized data from the scater object - the exprs slot contains log2(norm + 1)
mc_dc_normalized_with_scater <- as.matrix(exprs(sce))

#reverse the log2 + 1
norm_nolog2 <- (2^(mc_dc_normalized_with_scater)) -1

#natural log + 1
norm_ln2 <- log(norm_nolog2+ 1)

#load your original data into seurat
mc_dc <- CreateSeuratObject(raw.data = mc_dc, min.cells = 3, min.genes = 300, project = "MC_DC", meta.data = mc_dc_meta)

#load your newly normalized (natural log + 1) into the seurat.raw@data slot
mc_dc@data <- norm_ln2

mc_dc <- FindVariableGenes(object = mc_dc, mean.function = ExpMean, dispersion.function = LogVMR, 
                           x.low.cutoff = 0.06, x.high.cutoff = 2, y.cutoff = 0.5)



length(x = mc_dc@var.genes)
#1344
#scale the data
mc_dc <- ScaleData(object = mc_dc, vars.to.regress = c("nUMI", "percent.mito"))

#Run PCA
mc_dc <- RunPCA(object = mc_dc, pc.genes = mc_dc@var.genes, pcs.compute = 40, do.print = FALSE)
PCAPlot(object = mc_dc, dim.1 = 1, dim.2 = 2)


VizPCA(object = mc_dc, pcs.use = 1:2)
PrintPCA(object = mc_dc, pcs.print = 1:5)

mc_dc <- ProjectPCA(object = mc_dc, do.print = FALSE)
pdf("mc_dc_heatmap")
PCHeatmap(object = mc_dc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()

pdf("mc_dc_heatmap_20")
PCHeatmap(object = mc_dc, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()

mc_dc <- JackStraw(object = mc_dc, num.pc = 25, num.replicate = 100, do.print = FALSE)

JackStrawPlot(object = mc_dc, PCs = 1:25)

pdf("mc_dc_elbowplot")
PCElbowPlot(object = mc_dc,num.pc = 40)
dev.off()

#Find Clusters
mc_dc_12 <- FindClusters(object = mc_dc,reduction.type = "pca",dims.use = 1:11,resolution = 1.2, force.recalc= TRUE, save.SNN = FALSE)
mc_dc_12 <- RunTSNE(object = mc_dc_12, dims.use = 1:11,do.fast = TRUE)
TSNEPlot(object = mc_dc_12, do.label = TRUE)

# Merge clusters by renaming the identity class
mc_dc_12_merged <- RenameIdent(mc_dc_12, old.ident.name = 6, new.ident.name = 4 )

pdf("mc_dc_12_merged_tSNE.pdf")
TSNEPlot(mc_dc_12_merged, do.label = T)
dev.off()

mc_dc_12_merged <- StashIdent(mc_dc_12_merged, save.name = "merged.ident")
# rename the identities based on Slava’s grouping 
# divide cluster 8 into 8a and 8b manually 

7a
select.cells_c7a <- TSNEPlot(object = mc_dc_12_merged, do.identify = TRUE)

head(select.cells_c7a)
length(select.cells_c7a)
[1] 45

mc_dc_12_merged2 <- SetIdent(object = mc_dc_12_merged, cells.use = select.cells_c7a, ident.use = "7a")

#save the new identities for mc_dc_12_merged2
mc_dc_12_merged2 <- StashIdent(mc_dc_12_merged2, save.name = "merged2.ident")
merged2.ids <- c(0,1, 2, 3, 4, 5, 7,"7a",8)
new.merged2.ids <- c(5,2,4,1,6,3,"8b","8a",7)
mc_dc_12_merged2@ident <- plyr::mapvalues(x = mc_dc_12_merged2@ident, from = merged2.ids, to = new.merged2.ids)

mc_dc_12_merged2 <- StashIdent(mc_dc_12_merged2, save.name = "new.ident")

select.cells_3a <- TSNEPlot(object = mc_dc_12_merged2, do.identify = TRUE)
length(select.cells_3a)
[1] 29
mc_dc_12_merged2_subC3 <- SetIdent(object = mc_dc_12_merged2, cells.use = select.cells_3a, ident.use = "3a")

select.cells_3b <- TSNEPlot(object = mc_dc_12_merged2_subC3, do.identify = TRUE)
length(select.cells_3b)
[1] 16
mc_dc_12_merged2_subC3 <- SetIdent(object = mc_dc_12_merged2_subC3, cells.use = select.cells_3b, ident.use = "3b")

#save the new identities for mc_dc_12_merged2_subC3
mc_dc_12_merged2_subC3 <- StashIdent(mc_dc_12_merged2_subC3, save.name = "subC3.ident")

# Final tSNE for figures 1
pdf("mc_dc_12_merged2_sub3C_tSNE.pdf")
TSNEPlot(object = mc_dc_12_merged2_subC3, do.label = TRUE)
dev.off()

# Feature plots for figure 1
# this code will bring the red dots forward 
plot1 <- FeaturePlot(mc_dc_12_merged2_subC3, c("Adgre1", "C1qa", "Mafb", "Fcgr1"), cols.use = c("grey","red"), no.legend = TRUE)
pt.size = 1, pch.use = 16, no.legend = FALSE, do.return = TRUE, nCol = 2)
plot1 <- lapply(X = plot, function(p) { p$data <- p$data[order(p$data$gene),]; p})
cowplot::plot_grid(plotlist = plot1, ncol = ceiling(sqrt(length(plot1))))


plot2 <- FeaturePlot(mc_dc_12_merged2_subC3, c("Flt3", "Dpp4", "Xcr1", "Cd209a"), cols.use = c("grey","red"), no.legend = TRUE)
pt.size = 1, pch.use = 16, no.legend = FALSE, do.return = TRUE, nCol = 2)
plot2 <- lapply(X = plot2, function(p) { p$data <- p$data[order(p$data$gene),]; p})
cowplot::plot_grid(plotlist = plot2, ncol = ceiling(sqrt(length(plot2))))


plot3 <- FeaturePlot(mc_dc_12_merged2_subC3, c("Cdk1", "MKi67"), cols.use = c("grey","red"), no.legend = TRUE)
pt.size = 1, pch.use = 16, no.legend = FALSE, do.return = TRUE, nCol = 2)
plot3 <- lapply(X = plot, function(p) { p$data <- p$data[order(p$data$gene),]; p})
cowplot::plot_grid(plotlist = plot3, ncol = ceiling(sqrt(length(plot3))))


# Heatmap for supplementary figures)

# Change cluster numbers to re-order them for the paper 
mc_dc_12_merged2_subC3 <- SetAllIdent(mc_dc_12_merged2_subC3, id = "SubC3.ident")

current.cluster.ids <- c(1, 2,3,"3a","3b", 4, 5, 6, 7, "8a", "8b")
new.cluster.ids <- c(1,2,3,5,4,6,7,8,9,"10a", "10b")

mc_dc_12_merged2_subC3@ident <- plyr::mapvalues(x = mc_dc_12_merged2_subC3@ident, from = current.cluster.ids, to = new.cluster.ids)

mc_dc_12_merged2_subC3 <- StashIdent(mc_dc_12_merged2_subC3, save.name = "SubC2.ordered")

mc_dc_12_merged2_subC3 <- SetAllIdent(mc_dc_12_merged2_subC3, id = "SubC2.ordered")


TSNEPlot(object = mc_dc_12_merged2_subC3, do.label = TRUE, label.size =4)


#Find all markers
mc_dc_12_merged2_subC3.markers <- FindAllMarkers(mc_dc_12_merged2_subC3, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.2)


mc_dc_12_merged2_subC3.markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> T2.top20
mc_dc_12_merged2_subC3.markers %>% group_by(cluster) %>% top_n(65, avg_logFC) -> T2.top65

grep("Timd4", (T2.top65[which(T2.top65$cluster == 1), ])$gene)

pdf("mc_dc_12_merged2_subC3.top65markers_heatmap.pdf")
DoHeatmap(object = SubsetData(mc_dc_12_merged2_subC3, max.cells.per.ident = 100), genes.use = T2.top65$gene, cex.row = 3, slim.col.label = TRUE, remove.key = TRUE, group.cex = 12, group.spacing = 0.1)
dev.off()

pdf("mc_dc_12_merged2_subC3.top20markers_heatmap.pdf")
DoHeatmap(object = SubsetData(mc_dc_12_merged2_subC3, max.cells.per.ident = 100), genes.use = T2.top20$gene, cex.row = 2, slim.col.label = TRUE, remove.key = TRUE, group.cex = 12, group.spacing = 0.1)
dev.off()


#Violin Plots (for paper figure 1)

pdf("Figure1_Revesion2_vlnplots1.pdf")
VlnPlot(mc_dc_12_merged2_subC3, features.plot = c("Adgre1", "C1qa", "Mafb", "Fcgr1", "H2-Eb1", "Timd4", "Lyve1", "Folr2", "Igf1", "Ccr2", "Ace", "Plac8"), x.lab.rot = TRUE, size.x.use = 0, size.y.use = 7, size.title.use = 10, nCol = 3, same.y.lims = TRUE, y.log = TRUE, point.size.use= 0)
dev.off()

pdf("Figure1_Revesion2_vlnplots2.pdf")
VlnPlot(mc_dc_12_merged2_subC3, features.plot = c("Flt3", "Dpp4", "Xcr1", "Cd209a", "Cdk1", "Mki67", "Cdca3", "Birc5", "Il1b", "Ifit3", "Irf7", "Mmp12"), x.lab.rot = TRUE, size.x.use = 0, size.y.use = 7, size.title.use = 10, nCol = 3, same.y.lims = TRUE, y.log = TRUE, point.size.use= 0)
dev.off()



# Subset MCs and DC separately 

# (The section bellow was done based on the older cluster numbering… the cluster number change is above, )


dc.subset <- SubsetData(object = mc_dc_12_merged, ident.use =c(2,0,4))
mc.subset <- SubsetData(object = mc_dc_12_merged, ident.use =c(1,3,5))



dim(dc.subset@data)
[1] 12916   859

dim(mc.subset@data)
[1] 12916   641

# to perform additional rounds of clustering after SubsetData we recommend re-running FindVariableGenes() and ScaleData()

dc.subset <- FindVariableGenes(object = dc.subset, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = 0.5)

length(x = dc.subset@var.genes)
[1] 1468

mc.subset <- FindVariableGenes(object = mc.subset, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = 0.5)

length(x = mc.subset@var.genes)
[1] 1957

dc.subset <- ScaleData(object = dc.subset, vars.to.regress = c("nUMI", "percent.mito"))

mc.subset <- ScaleData(object = mc.subset, vars.to.regress = c("nUMI", "percent.mito"))



#Run PCA
dc.subset <- RunPCA(object = dc.subset, pc.genes = dc.subset@var.genes, pcs.compute = 40, do.print = FALSE)
mc.subset <- RunPCA(object = mc.subset, pc.genes = mc.subset@var.genes, pcs.compute =40, do.print = FALSE)


PCAPlot(object = dc.subset, dim.1 = 1, dim.2 = 2)
PCAPlot(object = mc.subset, dim.1 = 1, dim.2 = 2)

mc.subset <- ProjectPCA(object = mc.subset, do.print = FALSE)
dc.subset <- ProjectPCA(object = dc.subset, do.print = FALSE)


PCHeatmap(object = dc.subset, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = mc.subset, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)



PCHeatmap(object = dc.subset, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
#4

PCHeatmap(object = mc.subset, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
#6


PCElbowPlot(object = dc.subset,num.pc = 25)
# 5 PCs


PCElbowPlot(object = mc.subset,num.pc = 25)
# 9 PCs


# from the above methods, I think we will chose the first 6 significant 
#PCs for subsequent analysis of the mc.subset and the first 4 pcs for the dc.subset. 

dc.subset_res6 <- FindClusters(object = dc.subset,reduction.type = "pca",dims.use = 1:4,resolution = 0.6,save.SNN = FALSE, force.recalc = TRUE)
mc.subset_res6 <- FindClusters(object = mc.subset,reduction.type = "pca",dims.use = 1:9,resolution = 0.6,save.SNN = FALSE, force.recalc = TRUE)

dc.subset_res6 <- RunTSNE(object = dc.subset_6, dims.use = 1:4, resolution = 0.6, do.fast = TRUE)
mc.subset_res6 <- RunTSNE(object = mc.subset_res6, dims.use = 1:9,resolution = 0.6, do.fast = TRUE)

pdf("mc_subset_6_TSNEPlot.pdf")
TSNEPlot(object = mc.subset_res6, do.label = TRUE)
dev.off()

pdf("dc_subset_6_TSNEPlot.pdf")
TSNEPlot(object = dc.subset_res6, do.label = TRUE)
dev.off()

# Subset out cluster 2 of mc_subset_res6 and reculster


C2.subset <- SubsetData(object = mc.subset_res6, ident.use = 2 )

dim(C2.subset@data)
[1] 12916   109


C2.subset <- FindVariableGenes(object = C2.subset, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = C2.subset@var.genes)
[1] 2685

C2.subset <- ScaleData(object = C2.subset, vars.to.regress = c("nUMI", "percent.mito"))

#Run PCA
C2.subset <- RunPCA(object = C2.subset, pc.genes = C2.subset@var.genes, pcs.compute = 40, do.print = FALSE)

PCAPlot(object = C2.subset, dim.1 = 1, dim.2 = 2)

C2.subset <- ProjectPCA(object = C2.subset, do.print = FALSE)


PCHeatmap(object = C2.subset, pc.use = 1, cells.use = 100, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = C2.subset, pc.use = 1:20, cells.use = 100, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = C2.subset,num.pc = 30)

C2.subset_res6 <- FindClusters(object = C2.subset, reduction.type = "pca",dims.use = 1:20, resolution = 0.6,save.SNN = FALSE, force.recalc = TRUE)
C2.subset_res8 <- FindClusters(object = C2.subset, reduction.type = "pca",dims.use = 1:20, resolution = 0.8, save.SNN = FALSE, force.recalc = TRUE)
C2.subset_res9 <- FindClusters(object = C2.subset, reduction.type = "pca",dims.use = 1:20, resolution = 0.9, save.SNN = FALSE, force.recalc = TRUE)

C2.subset_res8<- RunTSNE(object = C2.subset_res8, dims.use = 1:20, resolution = 0.8, do.fast = TRUE)
C2.subset_res9<- RunTSNE(object = C2.subset_res9, dims.use = 1:20, resolution = 0.9, do.fast = TRUE)
C2.subset_res9<- RunTSNE(object = C2.subset_res9, dims.use = 1:20, do.fast = TRUE, add.iter = 2000) # didn’t change anything

setwd ("/Users/snejat/Documents/Epelman_filtered_gene_bc_matrices/mc_dc_results/mc_subset_res6")
pdf("C2_subset_9_TSNEPlot.pdf")
TSNEPlot(object = C2.subset_res9, do.label = TRUE)
dev.off()

Cluster comparisons for C2.subset_res9

#C1 vs C0
C2_subset_res9_C1vsC0markers =FindMarkers(C2.subset_res9 ,1,0,test.use = "MAST")
write.csv(C2_subset_res9_C1vsC0markers, file = "C2_mc_subset_C1vsC0markers.csv", sep = " ", col.names = TRUE, row.names = TRUE, quote = FALSE)


#C2 vs C0
C2_subset_res9_C2vsC0markers =FindMarkers(C2.subset_res9 ,2,0,test.use = "MAST")
write.csv(C2_subset_res9_C2vsC0markers, file = "C2_mc_subset_C2vsC0markers.csv", sep = " ", col.names = TRUE, row.names = TRUE, quote = FALSE)


#C2 vs C1
C2_subset_res9_C2vsC1markers =FindMarkers(C2.subset_res9 ,2,1,test.use = "MAST")
write.csv(C2_subset_res9_C2vsC1markers, file = "C2_mc_subset_C2vsC1markers.csv", sep = " ", col.names = TRUE, row.names = TRUE, quote = FALSE)

Produce a figure that subsets C2 into 2, 2a, 2b while maintaining other clusters


select.cells_c2a <- TSNEPlot(object = mc.subset_res6, do.identify = TRUE)
head(select.cells_c2a)
length(select.cells_c2a)
[1] 25

mc.subset_res6_subC2 <- SetIdent(object = mc.subset_res6, cells.use = select.cells_c2a, ident.use = "2a")

select.cells_c2b <- TSNEPlot(object = mc.subset_res6_subC2, do.identify = TRUE)
head(select.cells_c2b)
length(select.cells_c2b)
[1] 16

mc.subset_res6_subC2 <- SetIdent(object = mc.subset_res6_subC2, cells.use = select.cells_c2b, ident.use = "2b")

#save the new identities for mc_dc_12_merged2
mc.subset_res6_subC2 <- StashIdent(mc.subset_res6_subC2, save.name = "SubC2.ident")

# tSNE for figure 2 of manuscript
pdf("mc.subset_res6_subC2_tSNE.pdf")
TSNEPlot(object = mc.subset_res6_subC2, do.label = TRUE)
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
# generate a heatmap for  the top 20 markers

mc.subset_res6_subC2.markers <- FindAllMarkers(mc.subset_res6_subC2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

mc.subset_res6_subC2.markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> mc.subset_res6_subC2.markers_top20

mc.subset_res6_subC2.markers %>% group_by(cluster) %>% top_n(40, avg_logFC) -> mc.subset_res6_subC2.markers_top40


# HEAT MAP for figure S1 of the paper
setwd ( "/Users/snejat/Documents/Epelman_filtered_gene_bc_matrices/mc_dc_results/mc_subset_res6" )
pdf("mc.subset_6_subC2_top40markers_heatmap.pdf")
DoHeatmap(SubsetData(mc.subset_res6_subC2, max.cells.per.ident = 100) , genes.use = mc.subset_res6_subC2.markers_top40$gene, cex.row = 4, cex.col = 6, group.cex = 10, slim.col.label = TRUE, remove.key = TRUE, group.spacing = 0.10, group.label.rot = TRUE)
dev.off()

# Cluster comparisons ( to be used for Monocle)

#C1 vs C0
mc_subset_res6_subC2_C1vsC0markers =FindMarkers(mc.subset_res6_subC2 ,1,0,test.use = "MAST")
write.csv(mc_subset_res6_subC2_C1vsC0markers, file = "mc_subset_res6_subC2_C1vsC0markers.csv", row.names = TRUE, quote = FALSE)

list1 <- subset( mc_subset_res6_subC2_C1vsC0markers, mc_subset_res6_subC2_C1vsC0markers$p_val_adj < 0.05)


#C1 vs C2
mc_subset_res6_subC2_C1vsC2markers =FindMarkers(mc.subset_res6_subC2 ,1,2,test.use = "MAST")
write.csv(mc_subset_res6_subC2_C1vsC2markers, file = "mc_subset_res6_subC2_C1vsC2markers.csv", row.names = TRUE, quote = FALSE)

list2 <- subset( mc_subset_res6_subC2_C1vsC2markers, mc_subset_res6_subC2_C1vsC2markers$p_val_adj < 0.05)


#C1 vs C2a
mc_subset_res6_subC2_C1vsC2amarkers =FindMarkers(mc.subset_res6_subC2 ,1,"2a",test.use = "MAST")
write.csv(mc_subset_res6_subC2_C1vsC2amarkers, file = "mc_subset_res6_subC2_C1vsC2amarkers.csv", row.names = TRUE, quote = FALSE)

list3 <- subset( mc_subset_res6_subC2_C1vsC2amarkers, mc_subset_res6_subC2_C1vsC2amarkers$p_val_adj < 0.05)


#C1 vs C2b
mc_subset_res6_subC2_C1vsC2bmarkers =FindMarkers(mc.subset_res6_subC2 ,1,"2b",test.use = "MAST")
write.csv(mc_subset_res6_subC2_C1vsC2bmarkers, file = "mc_subset_res6_subC2_C1vsC2bmarkers.csv", row.names = TRUE, quote = FALSE)

list4 <- subset( mc_subset_res6_subC2_C1vsC2bmarkers, mc_subset_res6_subC2_C1vsC2bmarkers$p_val_adj < 0.05)


#C2 vs C0
mc_subset_res6_subC2_C2vsC0markers =FindMarkers(mc.subset_res6_subC2 ,2,0,test.use = "MAST")
write.csv(mc_subset_res6_subC2_C2vsC0markers, file = "mc_subset_res6_subC2_C2vsC0markers.csv", row.names = TRUE, quote = FALSE)

list5 <- subset( mc_subset_res6_subC2_C2vsC0markers, mc_subset_res6_subC2_C2vsC0markers$p_val_adj < 0.05)


#C2 vs C2a
mc_subset_res6_subC2_C2vsC2amarkers =FindMarkers(mc.subset_res6_subC2 ,2,"2a",test.use = "MAST")
write.csv(mc_subset_res6_subC2_C2vsC2amarkers, file = "mc_subset_res6_subC2_C2vsC2amarkers.csv", row.names = TRUE, quote = FALSE)

list6 <- subset( mc_subset_res6_subC2_C2vsC2amarkers, mc_subset_res6_subC2_C2vsC2amarkers$p_val_adj < 0.05)



#C2 vs C2b
mc_subset_res6_subC2_C2vsC2bmarkers =FindMarkers(mc.subset_res6_subC2 ,2,"2b",test.use = "MAST")
write.csv(mc_subset_res6_subC2_C2vsC2bmarkers, file = "mc_subset_res6_subC2_C2vsC2bmarkers.csv", row.names = TRUE, quote = FALSE)

list7 <- subset( mc_subset_res6_subC2_C2vsC2bmarkers, mc_subset_res6_subC2_C2vsC2bmarkers$p_val_adj < 0.05)


#C2a vs C0
mc_subset_res6_subC2_C2avsC0markers =FindMarkers(mc.subset_res6_subC2 ,"2a",0,test.use = "MAST")
write.csv(mc_subset_res6_subC2_C2avsC0markers, file = "mc_subset_res6_subC2_C2avsC0markers.csv", row.names = TRUE, quote = FALSE)

list8 <- subset( mc_subset_res6_subC2_C2avsC0markers , mc_subset_res6_subC2_C2avsC0markers $p_val_adj < 0.05)


#C2a vs C2b
mc_subset_res6_subC2_C2avsC2bmarkers =FindMarkers(mc.subset_res6_subC2 ,"2a","2b",test.use = "MAST")
write.csv(mc_subset_res6_subC2_C2avsC2bmarkers, file = "mc_subset_res6_subC2_C2avsC2bmarkers.csv", row.names = TRUE, quote = FALSE)

list9 <- subset( mc_subset_res6_subC2_C2avsC2bmarkers , mc_subset_res6_subC2_C2avsC2bmarkers $p_val_adj < 0.05)

#C2b vs C0
mc_subset_res6_subC2_C2bvsC0markers =FindMarkers(mc.subset_res6_subC2 ,"2b",0,test.use = "MAST")
write.csv(mc_subset_res6_subC2_C2bvsC0markers, file = "mc_subset_res6_subC2_C2bvsC0markers.csv", row.names = TRUE, quote = FALSE)

list10 <- subset(mc_subset_res6_subC2_C2bvsC0markers , mc_subset_res6_subC2_C2bvsC0markers $p_val_adj < 0.05)


# compile a list of statistically significant gene list from each cluster comparison above 
#( both up regulated and down regulated). This is to be used in Monocle for Psudotime Trajectory. 

mc_subset_res6_subC2_genelist <- append(row.names(list1), row.names(list2))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list3))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list4))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list5))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list6))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list7))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list8))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list9))
mc_subset_res6_subC2_genelist <- append(mc_subset_res6_subC2_genelist, row.names(list10))

length(mc_subset_res6_subC2_genelist )
[1] 1113


List <-mc_subset_res6_subC2_genelist[!duplicated(mc_subset_res6_subC2_genelist)]
length(List)
[1] 537
write.csv(List , file = "List .csv", quote = FALSE)

saveRDS (mc.subset_res6_subC2, file = "mc.subset_res6_subC2.rds")

# Make another list that only includes C2, C2a and C2b
mc_subset_res6_subC2_C2_C2a_C2b_genelist <- append(row.names(list6), row.names(list7))
mc_subset_res6_subC2_C2_C2a_C2b_genelist <- append(mc_subset_res6_subC2_C2_C2a_C2b_genelist, row.names(list9))

length(mc_subset_res6_subC2_C2_C2a_C2b_genelist)
[1] 58


List_C2_C2a_C2b <-mc_subset_res6_subC2_C2_C2a_C2b_genelist[!duplicated(mc_subset_res6_subC2_C2_C2a_C2b_genelist)]
length(List)
[1] 53

List_C2_C2a_C2b <- as.data.frame(List_C2_C2a_C2b)
write.csv(List_C2_C2a_C2b , file = "List_C2_C2a_C2b .csv", quote = FALSE)

# Subset out Clusters 2, 2a and 2b and do a Monocle trajectory on those 
dim(mc.subset_res6_subC2@data)
[1] 12916   641

mc.subset_res6_subC2_C2_C2a_C2b <- SubsetData(mc.subset_res6_subC2, ident.remove = c(0,1), subset.raw = TRUE)

dim(mc.subset_res6_subC2_C2_C2a_C2b@data)
[1] 12916   109

saveRDS (mc.subset_res6_subC2_C2_C2a_C2b, file = "mc.subset_res6_subC2_C2_C2a_C2b.rds")

############################################################
# Monocel Trajectory analysis
mc.subset_res6_subC2 <- readRDS("mc.subset_res6_subC2.rds")
packageVersion('monocle')
# ‘2.8.0’

dim(mc.subset_res6_subC2@raw.data)
#12744  1590

dim(mc.subset_res6_subC2@data)
# 12916   641

data6<-mc.subset_res6_subC2@raw.data[ , colnames(mc.subset_res6_subC2@data)]
dim(data6)
# 12744   641

#pheno data
cell.ids<-colnames(data6)
clust6<-as.data.frame(mc.subset_res6_subC2@ident)
cell.information <-cbind(cell.ids, clust6)

colnames(cell.information) <- c("cell.ids", "clust6")

pd6 <-new("AnnotatedDataFrame", data = data.frame(cell.information))

rownames(pd6)<-pd6$cell.ids

#feature data
gene.annot <-as.character(rownames(data6))
gene_short_name<-as.character(rownames(data6))
gene.annot2<-cbind(gene.annot, gene_short_name)

fd6<-new("AnnotatedDataFrame", data = data.frame(gene.annot2))
rownames(fd6)<-fd6$gene.annot


# create the DE gene list based on which the trajectory is formed
# FindMarkers for 2 clusters agaist each other at a time, then created genelist by appending the row 
# names of all the cluster comparison files, then removed duplicates (above).
# the file for this gene list is called List.csv

List <- read.csv("List.csv")
head(List)
#X     x
#1 1 Folr2
#2 2 F13a1
#3 3  Cbr2
#4 4   Pf4
#5 5 Lyve1
#6 6 Sepp1

ordering_genes <- List[,2]
head(ordering_genes)

mc_6_repeat <- setOrderingFilter(mc_6_repeat, ordering_genes )

mc_6_repeat <- estimateSizeFactors(mc_6_repeat)

mc_6_repeat <- reduceDimension(mc_6_repeat, max_components = 2,num_dim = 5,
                               method = 'DDRTree')

mc_6_repeat <- orderCells(mc_6_repeat)

pdf("mc_6_repeat_Pseudotime_bycluster.pdf")
plot_cell_trajectory(mc_6_repeat, color_by = "clust6")
dev.off()

######### Just to try the "import" function in Monocle ############
# did the same as above but use “import” for mc_subset_res6_subC2
#With a gene list created from cluster comparisons of mc_subset_res6_subC2 (same gene list as above)

mc_6 <- importCDS(mc.subset_res6_subC2, import_all = TRUE)

ordering_genes <- List[,2]
head(ordering_genes)

mc_6<- setOrderingFilter(mc_6, ordering_genes )


mc_6 <- estimateSizeFactors(mc_6)
mc_6 <- estimateDispersions(mc_6, relative_expr = TRUE, min_cells_detected =1, remove_outliers = TRUE)

mc_6 <- reduceDimension(mc_6, max_components = 2, num_dim = 5,
                        method = 'DDRTree')


mc_6 <- orderCells(mc_6)

pdf("mc.subset_res6_subC2_Pseudotime_bycluster_withgenelist.pdf")
plot_cell_trajectory(mc_6, color_by = "SubC2.ident")
dev.off()

pdf("mc.subset_res6_subC2_Pseudotime_byState_withgenelist.pdf")
plot_cell_trajectory(mc_6, color_by = "State")
dev.off()

#  ***  either way, results were the same ***  #

# to set a root state

GM_state <- function(mc_6){
  if (length(unique(pData(mc_6)$State)) > 1){
    T0_counts <- table(pData(mc_6)$State, pData(mc_6)$SubC2.ident)[,"2b"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}


mc_6 <- orderCells(mc_6, root_state = GM_state(mc_6))

pdf("mc.subset_res6_subC2_Pseudotime_bypsudotime_withgenelist.pdf")
plot_cell_trajectory(mc_6, color_by = "Pseudotime")
dev.off()

plot_cell_trajectory(mc_6, color_by = "SubC2.ident", show_tree = TRUE , markers = "Timd4", use_color_gradient = TRUE)

# plot each branch separately 
pdf("mc.subset_res6_subC2_Pseudotime_bycluster_withgenelist_separately.pdf")
plot_cell_trajectory(mc_6, color_by = "SubC2.ident") +
  facet_wrap(~SubC2.ident, nrow = 1)
dev.off()

Master_gene_list <- row.names(subset(fData(mc_6), gene_short_name %in% c("Retnla", "Lyve1", "Ccr2", "C1qa", "Plac8","Cx3cr1","Ly6c2","Ace","Ccl4", "Ccl9", "Ccl2", "Ccl3", "Ccl6", "CCl9", "Bcl2", "Itgax", "Runx3", "Igf1", "H2_Ab1", "H2_Aa", "Cd74", "Runx1", "PF4", "Spi1", "Csf1r", "C1qb", "Dnajb1", "Jund", "Anxa1", "Anxa2", "Klf2", "Klf4", "Fos", "Jun", "Junb", "Nfkbia", "Mrc1", "Cd169", "Tlr4", "Tlr2", "Lyz1", "Lyz2", "Sirpb1c", "H2-Eb1", "Il1b", "Mmp12","Cd9","Mmp13", "Maf", "Mafb", "Ifit1", "Ifit1bl1", "Ifit3", "Irf7" ,"S100a6", "Hspa1a", "Cxcl2", "S100a4", "Tmsb10", "Ccl12", "Ifitm3", "Apoe", "Folr2", "Lgals1", "Il1b","Sepp1","Wfdc17", "Timd4", "F13a1", "Sepp1", "Ccl24", "Gas6", "Ltc4s", "Cfh", "Fxyd2", "Ninj1", "Cd163", "Pltp", "Ifitm2", "Igfbp4", "Gngt2", "Itm2b", "Mrc1", "H2-DMb1", "Timp2", "Fcrls", "Fcgrt", "Maf", "Blvrb", "Pepd", "Fcna", "Cfp", "Psap", "Serpinb6a", "Rnase4", "Hexb", "Lgals1", "H2-DMa", "Pmp22", "Cd72", "Ctsd", "Smagp", "Snx2","Emp3", "C4b")))

pdf("Master_gene_list_heatmap.pdf")

plot_pseudotime_heatmap(mc_6[Master_gene_list,],
                        cluster_rows = TRUE,
                        cores = 1,
                        show_rownames = TRUE)
dev.off()

#Plot genes in Peseudotime
Genes_for_figure2 <- row.names(subset(fData(mc_6), gene_short_name %in% c("Plac8", "Ifit1","Irf7", "Cd9", "Ccr2", "Ly6c2","Klf4","Fos", "Jun", "Junb","Irf8", "Cx3cr1", "Runx3", "H2-DMa", "H2-DMb1", "H2-Eb1", "Spi1", "Igf1","Klf2", "Lyve1", "Timd4", "Folr2", "Retnla")))
plot_genes_in_pseudotime(mc_6[Genes_for_figure2,], color_by="SubC2.ident") + facet_wrap( ~ feature_label, scales= "free_y", ncol = 4 )

#######
# MONOCLE analysis of Clusters 2_2a_2b of the mc.subset_res6_subC2 data
# this object is a subset of the above object. 


mc_6_C2_C2a_C2b <- importCDS(mc.subset_res6_subC2_C2_C2a_C2b, import_all = TRUE)


# creat a new gene list (findmarker between the clusters of the mc.subset_res6_subC2_C2_C2a_C2b dataset)

List_C2_C2a_C2b <- read.csv("List_C2_C2a_C2b.csv")
head(List_C2_C2a_C2b)
List_C2_C2a_C2b
#1            Irf7
#2          Phf11b
#3          Ifit3b
#4            Zbp1
#5           Ifit3
#6            Bst2


ordering_genes <- List_C2_C2a_C2b[ ,1]
head(ordering_genes)

mc_6_C2_C2a_C2b<- setOrderingFilter(mc_6_C2_C2a_C2b, ordering_genes )


mc_6_C2_C2a_C2b <- estimateSizeFactors(mc_6_C2_C2a_C2b)

mc_6_C2_C2a_C2b <- reduceDimension(mc_6_C2_C2a_C2b, max_components = 2, num_dim = 3,
                                   method = 'DDRTree')


mc_6_C2_C2a_C2b <- orderCells(mc_6_C2_C2a_C2b)


pdf("mc.subset_res6_subC2_C2_C2a_C2b_Pseudotime_bycluster.pdf")
plot_cell_trajectory(mc_6_C2_C2a_C2b, color_by = "SubC2.ident", theta = 180, cell_size = 2, show_branch_points = FALSE)
dev.off()

# to set a root state

GM_state <- function(mc_6_C2_C2a_C2b){
  if (length(unique(pData(mc_6_C2_C2a_C2b)$State)) > 1){
    T0_counts <- table(pData(mc_6_C2_C2a_C2b)$State, pData(mc_6_C2_C2a_C2b)$SubC2.ident)[,"2b"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}


mc_6_C2_C2a_C2b <- orderCells(mc_6_C2_C2a_C2b, root_state = GM_state(mc_6_C2_C2a_C2b))


# plot each branch separately 

pdf("mc.subset_res6_subC2_Pseudotime_C2_C2a_C2b_bycluster_separately.pdf")
plot_cell_trajectory(mc_6_C2_C2a_C2b, color_by = "SubC2.ident") +
  facet_wrap(~SubC2.ident, nrow = 1)
dev.off()


### Mpath analysis #######
### Extract data for Mpath analysis ###

rds <- readRDS("mc.subset_res6_subC2.rds")

scaled_data <- rds@scale.data
write.table(scaled_data,"scaled_data.txt",sep="\t",col.names=NA)
write.table(rds@meta.data,"meta_data.txt",sep="\t",col.names=NA)
write.table(rds@meta.data[,c("SubC2.ident","cell.names")],"meta_data_cluster_celltype.txt",sep="\t",col.names=NA)
write.table(scaled_data[rds@var.genes,],"scaled_data_variable_genes.txt",sep="\t",col.names=NA)


library(Mpath)

distMethod = "euclidean"
rpkmFile="scaled_data_variable_genes.txt"
baseName="scaled_data_variable_genes_cluster_euclidean"

neighbor_network <- build_network(exprs = rpkmFile, iflog2 = FALSE,
                                  landmark_cluster = "meta_data_cluster.txt", distMethod = distMethod, baseName = baseName,textSize=10)

trimmed_net <- trim_net(neighbor_network,textSize=10,
                        baseName = baseName,
                        method = "mst")
