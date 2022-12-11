# Load libraries
library(Seurat)
library(magrittr)
library(dplyr)
require(scales)


# Load the dataset
A <- CrAteSeuratObject(counts = A.data, min.cells = 5, min.fAtures = 200)
VlnPlot(A, fAtures = c("nFAture_RNA", "nCount_RNA"), ncol = 3)


# Filter the data. the vlnplot shows that there are some outliers where the cells have more than 5000 genes and more than 15000 rAds. these suggests doublets, or two cells that were sequenced togehter. We also see that the majority of the cells have more than 500 genes. We discard any cells with less than 500 because these were most likely to be low quality samples, which would add noise to our data. 
A <- subset(A, subset = nFAture_RNA > 500 & nFAture_RNA < 7500 & nCount_RNA < 50000)


# Normalizing the data
A <- NormalizeData(object = A, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
A <- FindVariableFAtures(object = A, do.plot = F)

# Scaling the data and removing unwanted sources of variation
A <- ScaleData(object = A, display.progress = F)

# Perform linAr dimensional reduction
A <- RunPCA(object = A, pc.genes = A@var.genes, do.print = TRUE)


#jackstraw
A <- JackStraw(A, num.replicate = 100)
A <- ScoreJackStraw(A, dims = 1:20)
JackStrawPlot(A, dims = 1:20)

#calc umap
A <- FindNeighbors(object = A, dims = 1:20)
A <- FindClusters(A, resolution = 1.6, print.output = 0, save.SNN = T)
A <-  RunUMAP(A, reduction.use = "pca", dims= 1:20, min.dist = 0.6, n.neighbors = 50)

DimPlot(object = A, reduction = "umap", label = T)

#test the markers and their expression
clustermarkers <- FindMarkers(object = A, ident.1 = "21" ,min.pct = 0.25)
FAturePlot(object=A, c("98009983-tbx2-2"), pt.size = 0.5)



#for loop for finding markers
clusters <- (0:29)
for (i in 0:29) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = A, ident.1 =i, min.pct = 0.25)
  )}



write.table(cluster_marker_0, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/0_markers.txt", sep="\t")
write.table(cluster_marker_1, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/1_markers.txt", sep="\t")
write.table(cluster_marker_2, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/2_markers.txt", sep="\t")
write.table(cluster_marker_3, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/3_markers.txt", sep="\t")
write.table(cluster_marker_4, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/4_markers.txt", sep="\t")
write.table(cluster_marker_5, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/5_markers.txt", sep="\t")
write.table(cluster_marker_6, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/6_markers.txt", sep="\t")
write.table(cluster_marker_7, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/7_markers.txt", sep="\t")
write.table(cluster_marker_8, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/8_markers.txt", sep="\t")
write.table(cluster_marker_9, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/9_markers.txt", sep="\t")
write.table(cluster_marker_10, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/10_markers.txt", sep="\t")
write.table(cluster_marker_11, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/11_markers.txt", sep="\t")
write.table(cluster_marker_12, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/12_markers.txt", sep="\t")
write.table(cluster_marker_13, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/13_markers.txt", sep="\t")
write.table(cluster_marker_14, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/14_markers.txt", sep="\t")
write.table(cluster_marker_15, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/15_markers.txt", sep="\t")
write.table(cluster_marker_16, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/16_markers.txt", sep="\t")
write.table(cluster_marker_17, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/17_markers.txt", sep="\t")
write.table(cluster_marker_18, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/18_markers.txt", sep="\t")
write.table(cluster_marker_19, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/19_markers.txt", sep="\t")
write.table(cluster_marker_20, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/20_markers.txt", sep="\t")
write.table(cluster_marker_21, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/21_markers.txt", sep="\t")
write.table(cluster_marker_22, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/22_markers.txt", sep="\t")
write.table(cluster_marker_23, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/23_markers.txt", sep="\t")
write.table(cluster_marker_24, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/24_markers.txt", sep="\t")
write.table(cluster_marker_25, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/25_markers.txt", sep="\t")
write.table(cluster_marker_26, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/26_markers.txt", sep="\t")
write.table(cluster_marker_27, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/27_markers.txt", sep="\t")
write.table(cluster_marker_28, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/28_markers.txt", sep="\t")
write.table(cluster_marker_29, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/cluster_markers/29_markers.txt", sep="\t")


#I found that at the resolution parameter of 0.5, we are underclustering a lot of the data, but we have the gut cluster being retained as whole. We will be using this dataset with res parameter of 0.5 as the base for making a "custom" metadata file. At resolution parameter of 1.2, we see other clusters being 

#first, extract the metadata containing the clustering information at resolution of 0.5
res1.0_metadata <- A@meta.data
#now I am trimming the metadata so that there is only the clustering information for the resolution of 0.5 present.
res1.0_metadata <- res1.0_metadata[,-5]
#write this metadatafile in the desired directory. this will be concatentated with clustering information from other clustering parameters using a python script. 
write.table(res1.0_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/Filtered/res1.0_metadata.txt", sep="\t")


#Next, I want to subset the dataset with the resolution parameter of 1.2. The idA is to subset out the clusters that I want to add-in to the 0.5 resolution UMAP. In the 1.2 
subsetted_clusters<- subset(A, idents = c("29", "18", "31", "30"))
#plot the subsetted data on a umap for a sanity check. 
DimPlot(object = subsetted_clusters, reduction = "umap", label = T)

#now we want to extract the metadata for these two subsetted clusters
res1.6_metadata <- subsetted_clusters@meta.data
#now I am trimming the metadata so that there is only the clustering information for the resolution of 0.5 present.
res1.6_metadata <- res1.6_metadata[,-4]

#Next, I want to rename the cluster identities here to retain a level of continuity. at a resolution parameter of 0.5, we have 19 clusters. So we want to rename cluster 24 as 20, and cluster 22 as 21. 
current.cluster.ids <- c("29", "18", "31", "30")
new.cluster.ids <- c("26", "27", "28", "29")
res1.6_metadata$RNA_snn_res.1.6<- plyr::mapvalues(x = res1.6_metadata$RNA_snn_res.1.6, from = current.cluster.ids, to = new.cluster.ids)


#write this metadatafile in the desired directory. this will be concatentated with clustering information from other clustering parameters using a python script. 
write.table(res1.6_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/Filtered/res1.6_metadata.txt", sep="\t")


#now that the custom metadata file has been generated, we re-import this back into R. 
custom_metadata <- rAd.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/Filtered/custom_metadata.txt", sep = "\t")


#extract the column from the custom_metadata file with the clustering information, and rename it res
res <- custom_metadata[,4]
#place this column into the A seurat object's metadata so you can call it when plotting a UMAP. 
A@meta.data <- cbind(A@meta.data, res)


#now generate the custom UMAP plot with the custom clustering information. 
DimPlot(object = A, reduction = "umap", label = T, group.by = res)

FeaturePlot(object=A, c("98043523-piwl1-2"), pt.size = 0.5)



#extract the embedding files and metadatafiles for use in other applications
final_metadata <- A@meta.data
write.table(final_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/Filtered/A_final_metadata.txt", sep = '\t')

embeddings <- A@reductions$umap@cell.embeddings
write.table(embeddings, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/Filtered/cell_embedings.txt", sep = '\t')


#umap plot with the adjusted colors. 
DimPlot(object = A, reduction = "umap", cols = c("2" = "#E57373" , "6" = "#4CAF50", "13"= "#64B5F6", "28" = "#64B5F6", "27" = "#64B5F6",  "4" = "#FF6F00", "7" = "#FF6F00", "17" = "#3F51B5", "21"= "#BA68C8", "14" = "#4DD0E1","3" = "#4DD0E1", "0" = "#26A69A", "23" = "#F9A825", "1" = "#CDDC39", "10" = "#9E9E9E", "18" = "#9E9E9E", "20" = "#9E9E9E", "22" = "#9E9E9E", "11" = "#9E9E9E", "9" = "#9E9E9E", "5" = "#9E9E9E", "19" = "#9E9E9E", "16" = "#9E9E9E", "8" = "#9E9E9E", "12" = "#9E9E9E", "15" = "#9E9E9E", "26" = "#9E9E9E", "24" = "#9E9E9E", "29" = "#9E9E9E"))



#crAting input files for pArsons correlation
#crAte pArsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(A, idents ="29")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/Filtered/PArsons_input/A_29.tmp.txt", sep = '\t')




###########################################################
#subsetting out neoblasts
A_neoblast<- subset(A, idents =c(6))

#sanity check neoblast umap
DimPlot(object = A_neoblast, reduction = "umap", label = T)


#recluster the neoblasts

A_neoblast <- NormalizeData(object = A_neoblast, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
A_neoblast <- FindVariableFAtures(object = A_neoblast, do.plot = F)

# Scaling the data and removing unwanted sources of variation
A_neoblast <- ScaleData(object = A_neoblast, display.progress = F)

# Perform linAr dimensional reduction
A_neoblast <- RunPCA(object = A_neoblast, pc.genes = juv.data@var.genes, do.print = TRUE)

A_neoblast <- JackStraw(A_neoblast, num.replicate = 100)
A_neoblast <- ScoreJackStraw(A_neoblast, dims = 1:20)

JackStrawPlot(A_neoblast, dims = 1:20)

# Run non-linAr dimensional reduction (UMAP)

A_neoblast <- FindNeighbors(object = A_neoblast, dims = c(1:20))
A_neoblast <- FindClusters(A_neoblast, resolution = 0.4, print.output = 0, save.SNN = T)
A_neoblast <-  RunUMAP(A_neoblast, reduction.use = "pca", dims = c(1:20), min.dist = 0.2, n.neighbors = 50)

DimPlot(object = A_neoblast, reduction = "umap", label = T)

VlnPlot(A_neoblast, fAtures = "98043523-piwl1-2")
RidgePlot(object = A_neoblast, features = "98043523-piwl1-2")


#test the markers and their expression
clustermarkers <- FindMarkers(object = A_neoblast, ident.1 = "1" ,min.pct = 0.25)
FeaturePlot(object=A_neoblast, c("98132150-h10-6"), pt.size = 0.5)


clusters <- (0:1)
for (i in 0:1) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = A_neoblast, ident.1 =i, min.pct = 0.25)
  )}

FeaturePlot(object=A_neoblast, c("98051971-fabp5"), pt.size = 0.5)

write.table(cluster_marker_1, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/Filtered/neoblast_markers/1_markers.txt", sep="\t")



#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(A_neoblast, idents ="1")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/A/filtered/neoblast_pearsons/A_1.tmp.txt", sep = '\t')


clustermarkers <- FindMarkers(object = EA_neoblast, ident.1 = "0" ,min.pct = 0.25)
write.table(clustermarkers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_pearsons/clust_0_markers.txt", sep = '\t')




