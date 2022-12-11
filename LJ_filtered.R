# Load libraries
library(Seurat)
library(magrittr)
library(dplyr)
require(scales)


# Load the dataset

LJ <- CreateSeuratObject(counts = LJ_matrix, min.cells = 5, min.features = 200)
VlnPlot(LJ, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)


# Filter the data. the vlnplot shows that there are some outliers where the cells have more than 5000 genes and more than 15000 reads. these suggests doublets, or two cells that were sequenced togehter. We also see that the majority of the cells have more than 500 genes. We discard any cells with less than 500 because these were most likely to be low quality samples, which would add noise to our data. 
LJ <- subset(LJ, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)


# Normalizing the data
LJ <- NormalizeData(object = LJ, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
LJ <- FindVariableFeatures(object = LJ, do.plot = F)

# Scaling the data and removing unwanted sources of variation
LJ <- ScaleData(object = LJ, display.progress = F)

# Perform linear dimensional reduction
LJ <- RunPCA(object = LJ, pc.genes = LJ@var.genes, do.print = TRUE)


#jackstraw
LJ <- JackStraw(LJ, num.replicate = 100)
LJ <- ScoreJackStraw(LJ, dims = 1:20)
JackStrawPlot(LJ, dims = 1:20)

#calc umap
LJ <- FindNeighbors(object = LJ, dims = 1:20)
LJ <- FindClusters(LJ, resolution = 0.5, print.output = 0, save.SNN = T)
LJ <-  RunUMAP(LJ, reduction.use = "pca", dims= 1:20, min.dist = 0.6, n.neighbors = 50)

DimPlot(object = LJ, reduction = "umap", label = T)


#test the markers and their expression
clusters <- (0:19)
for (i in 0:19) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = LJ, ident.1 =i, min.pct = 0.25)
  )}

write.table(cluster_marker_0, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/0_markers.txt", sep="\t")
write.table(cluster_marker_1, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/1_markers.txt", sep="\t")
write.table(cluster_marker_2, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/2_markers.txt", sep="\t")
write.table(cluster_marker_3, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/3_markers.txt", sep="\t")
write.table(cluster_marker_4, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/4_markers.txt", sep="\t")
write.table(cluster_marker_5, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/5_markers.txt", sep="\t")
write.table(cluster_marker_6, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/6_markers.txt", sep="\t")
write.table(cluster_marker_7, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/7_markers.txt", sep="\t")
write.table(cluster_marker_8, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/8_markers.txt", sep="\t")
write.table(cluster_marker_9, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/9_markers.txt", sep="\t")
write.table(cluster_marker_10, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/10_markers.txt", sep="\t")
write.table(cluster_marker_11, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/11_markers.txt", sep="\t")
write.table(cluster_marker_12, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/12_markers.txt", sep="\t")
write.table(cluster_marker_13, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/13_markers.txt", sep="\t")
write.table(cluster_marker_14, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/14_markers.txt", sep="\t")
write.table(cluster_marker_15, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/15_markers.txt", sep="\t")
write.table(cluster_marker_16, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/16_markers.txt", sep="\t")
write.table(cluster_marker_17, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/17_markers.txt", sep="\t")
write.table(cluster_marker_18, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/18_markers.txt", sep="\t")
write.table(cluster_marker_19, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_cluster_markers/19_markers.txt", sep="\t")






FeaturePlot(object=LJ, c("98022905-h2ax-2"), pt.size = 0.5)



#umap plot with the adjusted colors. 
DimPlot(object = LJ, reduction = "umap", cols = c("4" = "#E57373" , "1" = "#4CAF50", "8"= "#64B5F6","9" = "#64B5F6", "12" = "#64B5F6", "0" = "#FF6F00", "2" = "#3F51B5", "13" = "#3F51B5", "7"= "#BA68C8", "6" = "#4DD0E1", "18" = "#26A69A", "15" = "#BA68C8", "5" = "#F9A825", "3" = "#CDDC39", "19" = "#9E9E9E", "11" = "#9E9E9E", "17" = "#9E9E9E", "10" = "#9E9E9E", "14" = "#9E9E9E", "16" = "#9E9E9E"))


LJ@meta.data <- LJ@meta.data[,-5]


#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(LJ, idents ="19")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/pearsons_input/LJ_19.tmp.txt", sep = '\t')



#make sure the active idents are from res in the neoblast subset. 
#LJ <- SetIdent(LJ, value = LJ@meta.data$res)
LJ@active.ident


###########################################################
#subsetting out neoblasts
LJ_neoblast<- subset(LJ, idents ="1")

#sanity check neoblast umap
DimPlot(object = LJ_neoblast, reduction = "umap", label = T)


#recluster the neoblasts

LJ_neoblast <- NormalizeData(object = LJ_neoblast, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
LJ_neoblast <- FindVariableFeatures(object = LJ_neoblast, do.plot = F)

# Scaling the data and removing unwanted sources of variation
LJ_neoblast <- ScaleData(object = LJ_neoblast, display.progress = F)

# Perform linear dimensional reduction
LJ_neoblast <- RunPCA(object = LJ_neoblast, pc.genes = juv.data@var.genes, do.print = TRUE)

LJ_neoblast <- JackStraw(LJ_neoblast, num.replicate = 100)
LJ_neoblast <- ScoreJackStraw(LJ_neoblast, dims = 1:20)

JackStrawPlot(LJ_neoblast, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

LJ_neoblast <- FindNeighbors(object = LJ_neoblast, dims = c(1:13,15:17))
LJ_neoblast <- FindClusters(LJ_neoblast, resolution = 1.0, print.output = 0, save.SNN = T)
LJ_neoblast <-  RunUMAP(LJ_neoblast, reduction.use = "pca", dims = c(1:13,15:17), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = LJ_neoblast, reduction = "umap", label = T)

VlnPlot(LJ_neoblast, features = "98043523-piwl1-2")
RidgePlot(object = LJ_neoblast, features = "98043523-piwl1-2")


#test the markers and their expression
clustermarkers <- FindMarkers(object = LJ_neoblast, ident.1 = "4" ,min.pct = 0.25)
FeaturePlot(object=LJ_neoblast, c("98043523-piwl1-2"), pt.size = 0.5)


clusters <- (0:7)
for (i in 0:7) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = LJ_neoblast, ident.1 =i, min.pct = 0.25)
  )}


RidgePlot(object = LJ_neoblast, features = "98043523-piwl1-2")


write.table(cluster_marker_5, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/2022_pearsons_neoblast/neoblast_cluster_markers/5_markers.txt", sep="\t")




#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(LJ_neoblast, idents ="6")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/2022_pearsons_neoblast/LJ_6.tmp.txt", sep = '\t')


clustermarkers <- FindMarkers(object = LJ_neoblast, ident.1 = "0" ,min.pct = 0.25)
write.table(clustermarkers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_pearsons/clust_0_markers.txt", sep = '\t')



write.table(LJ@assays$RNA@data, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/LJ_neoblast_matrix.txt", sep = '\t')


