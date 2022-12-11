library(Seurat)

x17dpa_matrix <- read.table("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/sc.counts_17dpa.matrix", sep = "\t")


reg_17dpa <- CreateSeuratObject(counts = x17dpa_matrix, min.cells = 5, min.features =200)


# Filter the data
####script for embedding the post embryonic datasets myself
# Filter the data
####script for embedding the post embryonic datasets myself
VlnPlot(reg_17dpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
reg_17dpa <- subset(reg_17dpa, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000)



# Normalizing the data
reg_17dpa <- NormalizeData(object = reg_17dpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_17dpa <- FindVariableFeatures(object =reg_17dpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_17dpa <- ScaleData(object = reg_17dpa, display.progress = F)

# Perform linear dimensional reduction
reg_17dpa <- RunPCA(object =reg_17dpa, pc.genes = juv.data@var.genes, do.print = TRUE)

reg_17dpa <- JackStraw(reg_17dpa, num.replicate = 100)
reg_17dpa <- ScoreJackStraw(reg_17dpa, dims = 1:20)

JackStrawPlot(reg_17dpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

reg_17dpa <- FindNeighbors(object = reg_17dpa, dims = c(1:20))
reg_17dpa <- FindClusters(reg_17dpa, resolution = 1.5, print.output = 0, save.SNN = T)
reg_17dpa <-  RunUMAP(reg_17dpa, reduction.use = "pca", dims= c(1:20), min.dist = 0.5, n.neighbors = 50)

DimPlot(object = reg_17dpa, reduction = "umap",label = T)


cluster0.markers <- FindMarkers(object = reg_17dpa, ident.1 =0, min.pct = 0.25)


FeaturePlot(reg_17dpa, c("98048987-catl2-4"), pt.size = 0.5)



VlnPlot(reg_17dpa, features = "98043523-piwl1-2")
RidgePlot(object =reg_17dpa, features = "98043523-piwl1-2")




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
neoblast_17dpa<- subset(reg_17dpa, idents =3)
DimPlot(object = neoblast_17dpa, reduction = "umap",label = T)

neoblast_17dpa <- NormalizeData(object = neoblast_17dpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
neoblast_17dpa <- FindVariableFeatures(object = neoblast_17dpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
neoblast_17dpa <- ScaleData(object = neoblast_17dpa, display.progress = F)

# Perform linear dimensional reduction
neoblast_17dpa <- RunPCA(object = neoblast_17dpa, pc.genes = juv.data@var.genes, do.print = TRUE)

neoblast_17dpa <- JackStraw(neoblast_17dpa, num.replicate = 100)
neoblast_17dpa <- ScoreJackStraw(neoblast_17dpa, dims = 1:20)

JackStrawPlot(neoblast_17dpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

neoblast_17dpa <- FindNeighbors(object = neoblast_17dpa, dims = c(1:5,7:9,13))
neoblast_17dpa <- FindClusters((neoblast_17dpa), resolution = 1.0, print.output = 0, save.SNN = T)
neoblast_17dpa <-  RunUMAP(neoblast_17dpa, reduction.use = "pca", dims= c(1:5,7:9,13), min.dist = 0.2, n.neighbors = 50)

DimPlot(object = reg_17dpa, reduction = "umap", cols = c("26" = "#E57373", "1"= "#E57373" ,"7"= "#E57373", "3" = "#4CAF50", "12"= "#64B5F6", "8" = "#64B5F6", "16" = "#64B5F6", "15" = "#64B5F6", "13" = "#FF6F00","11" = "#FF6F00", "21" = "#3F51B5", "4" = "#3F51B5","18"= "#BA68C8", "6" = "#4DD0E1","5" = "#4DD0E1", "9" = "#26A69A", "20" = "#BA68C8", "19" = "#F9A825", "0" = "#CDDC39","2" = "#CDDC39","24" = "#CDDC39", "14" = "#9E9E9E", "25" = "#9E9E9E", "27" = "#9E9E9E", "22" = "#9E9E9E", "23" = "#9E9E9E", "17" = "#9E9E9E", "10" = "#9E9E9E"))

VlnPlot(neoblast_17dpa, features = "98043523-piwl1-2")
RidgePlot(object = neoblast_17dpa, features = "98043523-piwl1-2")


#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_17dpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_17dpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_17dpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_17dpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_17dpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_17dpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object =neoblast_17dpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object =neoblast_17dpa, ident.1 =7, min.pct = 0.25)
cluster8.markers <- FindMarkers(object =neoblast_17dpa, ident.1 =8, min.pct = 0.25)

FeaturePlot(neoblast_17dpa, c("98009983-tbx2-2"), pt.size = 0.5)


#writing the cluster markers into files
cluster0.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "0", min.pct = 0.25)
write.table(cluster0.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust0_markers.txt", sep="\t")

cluster1.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "1", min.pct = 0.25)
write.table(cluster1.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust1_markers.txt", sep="\t")


cluster2.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "2", min.pct = 0.25)
write.table(cluster2.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust2_markers.txt", sep="\t")


cluster3.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "3", min.pct = 0.25)
write.table(cluster3.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust3_markers.txt", sep="\t")


cluster4.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "4", min.pct = 0.25)
write.table(cluster4.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust4_markers.txt", sep="\t")


cluster5.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "5", min.pct = 0.25)
write.table(cluster5.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust5_markers.txt", sep="\t")


cluster6.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "6", min.pct = 0.25)
write.table(cluster6.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust6_markers.txt", sep="\t")


cluster7.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "7", min.pct = 0.25)
write.table(cluster7.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust7_markers.txt", sep="\t")

cluster8.markers <- FindMarkers(object = neoblast_17dpa, ident.1 = "8", min.pct = 0.25)
write.table(cluster8.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/17dpa/cluster4_subset/cluster_markers/clust8_markers.txt", sep="\t")



#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(neoblast_17dpa, idents ="6")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_reg_neoblast/pearsons/x17dpa_6N.tmp.txt", sep = '\t')

