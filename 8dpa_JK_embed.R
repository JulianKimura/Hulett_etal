library(Seurat)

x8dpa_matrix <- read.table("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/sc.counts_8dpa.matrix", sep = "\t")


reg_8dpa <- CreateSeuratObject(counts = x8dpa_matrix, min.cells = 5, min.features =200)

VlnPlot(reg_8dpa, features = c("nFeature_RNA", "nCount_RNA"))


# Filter the data
####script for embedding the post embryonic datasets myself

# Filter the data
####script for embedding the post embryonic datasets myself
VlnPlot(reg_8dpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
reg_8dpa <- subset(reg_8dpa, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000)


# Normalizing the data
reg_8dpa <- NormalizeData(object = reg_8dpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_8dpa <- FindVariableFeatures(object =reg_8dpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_8dpa <- ScaleData(object = reg_8dpa, display.progress = F)

# Perform linear dimensional reduction
reg_8dpa <- RunPCA(object =reg_8dpa, pc.genes = juv.data@var.genes, do.print = TRUE)

reg_8dpa <- JackStraw(reg_8dpa, num.replicate = 100)
reg_8dpa <- ScoreJackStraw(reg_8dpa, dims = 1:20)

JackStrawPlot(reg_8dpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

reg_8dpa <- FindNeighbors(object = reg_8dpa, dims = c(1:20))
reg_8dpa <- FindClusters(reg_8dpa, resolution = 1.5, print.output = 0, save.SNN = T)
reg_8dpa <-  RunUMAP(reg_8dpa, reduction.use = "pca", dims= c(1:20), min.dist = 0.5, n.neighbors = 50)

DimPlot(object = reg_8dpa, reduction = "umap",label = T)


cluster0.markers <- FindMarkers(object = reg_8dpa, ident.1 =0, min.pct = 0.25)

FeaturePlot(reg_8dpa, c("98048987-catl2-4"), pt.size = 0.5)



VlnPlot(reg_8dpa, features = "98043523-piwl1-2")
RidgePlot(object =reg_8dpa, features = "98043523-piwl1-2")




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
neoblast_8dpa<- subset(reg_8dpa, idents =2)
DimPlot(object = neoblast_8dpa, reduction = "umap",label = T)

neoblast_8dpa <- NormalizeData(object = neoblast_8dpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
neoblast_8dpa <- FindVariableFeatures(object = neoblast_8dpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
neoblast_8dpa <- ScaleData(object = neoblast_8dpa, display.progress = F)

# Perform linear dimensional reduction
neoblast_8dpa <- RunPCA(object = neoblast_8dpa, pc.genes = juv.data@var.genes, do.print = TRUE)

neoblast_8dpa <- JackStraw(neoblast_8dpa, num.replicate = 100)
neoblast_8dpa <- ScoreJackStraw(neoblast_8dpa, dims = 1:20)

JackStrawPlot(neoblast_8dpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

neoblast_8dpa <- FindNeighbors(object = neoblast_8dpa, dims = c(1:11,13,14,15,16))
neoblast_8dpa <- FindClusters((neoblast_8dpa), resolution = 0.7, print.output = 0, save.SNN = T)
neoblast_8dpa <-  RunUMAP(neoblast_8dpa, reduction.use = "pca", dims= c(1:11,13,14,15,16), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = reg_8dpa, reduction = "umap", cols = c("21" = "#E57373", "1"= "#E57373" , "2" = "#4CAF50", "6"= "#64B5F6", "23" = "#64B5F6", "13" = "#64B5F6", "8" = "#64B5F6", "11" = "#FF6F00", "19" = "#FF6F00", "7" = "#3F51B5", "22"= "#BA68C8", "12" = "#4DD0E1","0" = "#4DD0E1", "4" = "#4DD0E1", "16" = "#26A69A", "10" = "#F9A825", "25" = "#CDDC39","3" = "#CDDC39","14" = "#CDDC39", "9" = "#9E9E9E", "18" = "#9E9E9E", "17" = "#9E9E9E", "5" = "#9E9E9E", "27" = "#9E9E9E", "20" = "#9E9E9E", "26" = "#9E9E9E", "24" = "#9E9E9E", "15" = "#9E9E9E"))

VlnPlot(neoblast_8dpa, features = c("nFeature_RNA", "nCount_RNA"))


VlnPlot(neoblast_8dpa, features = "98043523-piwl1-2")
RidgePlot(object = neoblast_8dpa, features = "98043523-piwl1-2")


#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_8dpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_8dpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_8dpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_8dpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_8dpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_8dpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object =neoblast_8dpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object =neoblast_8dpa, ident.1 =7, min.pct = 0.25)
cluster8.markers <- FindMarkers(object =neoblast_8dpa, ident.1 =8, min.pct = 0.25)

FeaturePlot(neoblast_8dpa, c("98050114-h33-3"), pt.size = 0.5)


#writing the cluster markers into files
cluster0.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "0", min.pct = 0.25)
write.table(cluster0.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust0_markers.txt", sep="\t")

cluster1.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "1", min.pct = 0.25)
write.table(cluster1.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust1_markers.txt", sep="\t")


cluster2.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "2", min.pct = 0.25)
write.table(cluster2.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust2_markers.txt", sep="\t")


cluster3.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "3", min.pct = 0.25)
write.table(cluster3.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust3_markers.txt", sep="\t")


cluster4.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "4", min.pct = 0.25)
write.table(cluster4.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust4_markers.txt", sep="\t")


cluster5.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "5", min.pct = 0.25)
write.table(cluster5.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust5_markers.txt", sep="\t")


cluster6.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "6", min.pct = 0.25)
write.table(cluster6.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust6_markers.txt", sep="\t")


cluster7.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "7", min.pct = 0.25)
write.table(cluster7.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust7_markers.txt", sep="\t")

cluster8.markers <- FindMarkers(object = neoblast_8dpa, ident.1 = "8", min.pct = 0.25)
write.table(cluster8.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/8dpa/cluster2_subset/cluster_markers/clust8_markers.txt", sep="\t")



#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(neoblast_8dpa, idents ="5")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_reg_neoblast/pearsons/x8dpa_5N.tmp.txt", sep = '\t')
