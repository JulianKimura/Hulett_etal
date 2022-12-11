library(Seurat)

x29dpa_matrix <- read.table("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/sc.counts_29dpa.matrix", sep = "\t")

rm(list=ls(pattern="neoblast"))



reg_29dpa <- CreateSeuratObject(counts = x29dpa_matrix, min.cells = 5, min.features =200)


# Filter the data
####script for embedding the post embryonic datasets myself
VlnPlot(reg_29dpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
reg_29dpa <- subset(reg_29dpa, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 20000)



# Normalizing the data
reg_29dpa <- NormalizeData(object = reg_29dpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_29dpa <- FindVariableFeatures(object =reg_29dpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_29dpa <- ScaleData(object = reg_29dpa, display.progress = F)

# Perform linear dimensional reduction
reg_29dpa <- RunPCA(object =reg_29dpa, pc.genes = reg_29dpa@var.genes, do.print = TRUE)

reg_29dpa <- JackStraw(reg_29dpa, num.replicate = 100)
reg_29dpa <- ScoreJackStraw(reg_29dpa, dims = 1:20)

JackStrawPlot(reg_29dpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

reg_29dpa <- FindNeighbors(object = reg_29dpa, dims = c(1:20))
reg_29dpa <- FindClusters(reg_29dpa, resolution = 1.7, print.output = 0, save.SNN = T)
reg_29dpa <-  RunUMAP(reg_29dpa, reduction.use = "pca", dims= c(1:20), min.dist = 0.2, n.neighbors = 50)


DimPlot(object = reg_29dpa, reduction = "umap", cols = c("7" = "#E57373", "8"= "#E57373" ,"26"= "#E57373", "5" = "#4CAF50", "10"= "#64B5F6", "13" = "#64B5F6", "15" = "#64B5F6", "30" = "#64B5F6","24" = "#64B5F6","12" = "#64B5F6", "8" = "#FF6F00","4" = "#FF6F00", "11" = "#3F51B5", "17"= "#BA68C8", "3" = "#4DD0E1","9" = "#4DD0E1","22" = "#4DD0E1", "6" = "#26A69A", "19" = "#F9A825", "2" = "#CDDC39","1" = "#CDDC39", "28" = "#9E9E9E", "16" = "#9E9E9E", "23" = "#9E9E9E", "21" = "#9E9E9E", "0" = "#9E9E9E", "27" = "#9E9E9E", "25" = "#9E9E9E", "20" = "#9E9E9E", "29" = "#9E9E9E", "14" = "#9E9E9E"))


cluster15.markers <- FindMarkers(object = reg_29dpa, ident.1 =15, min.pct = 0.25)

FeaturePlot(reg_29dpa, c("98048987-catl2-4"), pt.size = 0.5)



VlnPlot(reg_29dpa, features = "98043523-piwl1-2")
RidgePlot(object =reg_29dpa, features = "98043523-piwl1-2")




#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- SubsetData(LJ, ident.use ="3")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/URD/REH_paper/RH_For_JK_URD_Data2/JK_embedding/LJ/subset_cluster9/LJ_3.tmp.txt", sep = '\t')




########################## SUBSETTING BASED ON NEOBLAST CLUSTER (CLUSTER 0)
#first, subset based on cluster
neoblast_29dpa<- subset(reg_29dpa, idents =5)
DimPlot(object = neoblast_29dpa, reduction = "umap",label = T)

neoblast_29dpa <- NormalizeData(object = neoblast_29dpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
neoblast_29dpa <- FindVariableFeatures(object = neoblast_29dpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
neoblast_29dpa <- ScaleData(object = neoblast_29dpa, display.progress = F)

# Perform linear dimensional reduction
neoblast_29dpa <- RunPCA(object = neoblast_29dpa, pc.genes = neoblast_29dpa@var.genes, do.print = TRUE)

neoblast_29dpa <- JackStraw(neoblast_29dpa, num.replicate = 100)
neoblast_29dpa <- ScoreJackStraw(neoblast_29dpa, dims = 1:20)

JackStrawPlot(neoblast_29dpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

neoblast_29dpa <- FindNeighbors(object = neoblast_29dpa, dims = c(1:11,13))
neoblast_29dpa <- FindClusters((neoblast_29dpa), resolution = 0.5, print.output = 0, save.SNN = T)
neoblast_29dpa <-  RunUMAP(neoblast_29dpa, reduction.use = "pca", dims= c(1:11,13), min.dist = 0.2, n.neighbors = 50)

DimPlot(object = neoblast_29dpa, reduction = "umap",label = T)

VlnPlot(neoblast_29dpa, features = "98043523-piwl1-2")
RidgePlot(object = neoblast_29dpa, features = "98043523-piwl1-2")


#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_29dpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_29dpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_29dpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_29dpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_29dpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_29dpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object =neoblast_29dpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object =neoblast_29dpa, ident.1 =7, min.pct = 0.25)
cluster8.markers <- FindMarkers(object =neoblast_29dpa, ident.1 =8, min.pct = 0.25)
cluster9.markers <- FindMarkers(object =neoblast_29dpa, ident.1 =9, min.pct = 0.25)

FeaturePlot(neoblast_29dpa, c("98009983-tbx2-2"), pt.size = 0.5)

write.table(neoblast_29dpa@assays$RNA@data, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/filtered/29dpa_neoblast_subset.txt", sep = "\t")

#writing the cluster markers into files
cluster0.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "0", min.pct = 0.25)
write.table(cluster0.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust0_markers.txt", sep="\t")

cluster1.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "1", min.pct = 0.25)
write.table(cluster1.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust1_markers.txt", sep="\t")


cluster2.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "2", min.pct = 0.25)
write.table(cluster2.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust2_markers.txt", sep="\t")


cluster3.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "3", min.pct = 0.25)
write.table(cluster3.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust3_markers.txt", sep="\t")


cluster4.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "4", min.pct = 0.25)
write.table(cluster4.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust4_markers.txt", sep="\t")


cluster5.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "5", min.pct = 0.25)
write.table(cluster5.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust5_markers.txt", sep="\t")


cluster6.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "6", min.pct = 0.25)
write.table(cluster6.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust6_markers.txt", sep="\t")


cluster7.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "7", min.pct = 0.25)
write.table(cluster7.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust7_markers.txt", sep="\t")

cluster8.markers <- FindMarkers(object = neoblast_29dpa, ident.1 = "8", min.pct = 0.25)
write.table(cluster8.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/29dpa/cluster5_subset/cluster_markers/clust8_markers.txt", sep="\t")





#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(neoblast_29dpa, idents ="4")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_reg_neoblast/pearsons/x29dpa_4N.tmp.txt", sep = '\t')

