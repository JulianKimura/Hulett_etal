# Load libraries
library(Seurat)
library(magrittr)
library(dplyr)
require(scales)


# Load the dataset
EA <- CreateSeuratObject(counts = EA.data, min.cells = 5, min.features = 200)
VlnPlot(EA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)


# Filter the data. the vlnplot shows that there are some outliers where the cells have more than 5000 genes and more than 15000 reads. these suggests doublets, or two cells that were sequenced togehter. We also see that the majority of the cells have more than 500 genes. We discard any cells with less than 500 because these were most likely to be low quality samples, which would add noise to our data. 
EA <- subset(EA, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 50000)


# Normalizing the data
EA <- NormalizeData(object = EA, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
EA <- FindVariableFeatures(object = EA, do.plot = F)

# Scaling the data and removing unwanted sources of variation
EA <- ScaleData(object = EA, display.progress = F)

# Perform linear dimensional reduction
EA <- RunPCA(object = EA, pc.genes = EA@var.genes, do.print = TRUE)


#jackstraw
EA <- JackStraw(EA, num.replicate = 100)
EA <- ScoreJackStraw(EA, dims = 1:20)
JackStrawPlot(EA, dims = 1:20)

#calc umap
EA <- FindNeighbors(object = EA, dims = 1:20)
EA <- FindClusters(EA, resolution = 1.6, print.output = 0, save.SNN = T)
EA <-  RunUMAP(EA, reduction.use = "pca", dims= 1:20, min.dist = 0.6, n.neighbors = 50)

DimPlot(object = EA, reduction = "umap", label = T, group.by = 'res')

#test the markers and their expression
clustermarkers <- FindMarkers(object = EA, ident.1 = "21" ,min.pct = 0.25)
FeaturePlot(object=EA, c("98009022-calu"), pt.size = 0.5)


#I found that at the resolution parameter of 0.5, we are underclustering a lot of the data, but we have the gut cluster being retained as whole. We will be using this dataset with res parameter of 0.5 as the base for making a "custom" metadata file. At resolution parameter of 1.2, we see other clusters being 

#first, extract the metadata containing the clustering information at resolution of 0.5
res1.0_metadata <- EA@meta.data
#now I am trimming the metadata so that there is only the clustering information for the resolution of 0.5 present.
res1.0_metadata <- res1.0_metadata[,-5]
#write this metadatafile in the desired directory. this will be concatentated with clustering information from other clustering parameters using a python script. 
write.table(res1.0_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/Filtered/res1.0_metadata.txt", sep="\t")


#Next, I want to subset the dataset with the resolution parameter of 1.2. The idea is to subset out the clusters that I want to add-in to the 0.5 resolution UMAP. In the 1.2 
subsetted_clusters<- subset(EA, idents = c("29", "18", "31", "30"))
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
write.table(res1.6_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/Filtered/res1.6_metadata.txt", sep="\t")


#now that the custom metadata file has been generated, we re-import this back into R. 
custom_metadata <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/Filtered/custom_metadata.txt", sep = "\t")


#extract the column from the custom_metadata file with the clustering information, and rename it res
res <- custom_metadata[,4]
#place this column into the EA seurat object's metadata so you can call it when plotting a UMAP. 
EA@meta.data <- cbind(EA@meta.data, res)



EA <- SetIdent(EA, value = EA@meta.data$res)
EA@active.ident

#now generate the custom UMAP plot with the custom clustering information. 
clusters <- (26:29)
for (i in 26:29) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = EA, ident.1 =i, min.pct = 0.25)
  )}


write.table(cluster_marker_0, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/0_markers.txt", sep="\t")
write.table(cluster_marker_1, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/1_markers.txt", sep="\t")
write.table(cluster_marker_2, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/2_markers.txt", sep="\t")
write.table(cluster_marker_3, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/3_markers.txt", sep="\t")
write.table(cluster_marker_4, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/4_markers.txt", sep="\t")
write.table(cluster_marker_5, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/5_markers.txt", sep="\t")
write.table(cluster_marker_6, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/6_markers.txt", sep="\t")
write.table(cluster_marker_7, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/7_markers.txt", sep="\t")
write.table(cluster_marker_8, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/8_markers.txt", sep="\t")
write.table(cluster_marker_9, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/9_markers.txt", sep="\t")
write.table(cluster_marker_10, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/10_markers.txt", sep="\t")
write.table(cluster_marker_11, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/11_markers.txt", sep="\t")
write.table(cluster_marker_12, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/12_markers.txt", sep="\t")
write.table(cluster_marker_13, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/13_markers.txt", sep="\t")
write.table(cluster_marker_15, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/15_markers.txt", sep="\t")
write.table(cluster_marker_16, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/16_markers.txt", sep="\t")
write.table(cluster_marker_17, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/17_markers.txt", sep="\t")
write.table(cluster_marker_18, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/18_markers.txt", sep="\t")
write.table(cluster_marker_19, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/19_markers.txt", sep="\t")
write.table(cluster_marker_20, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/20_markers.txt", sep="\t")
write.table(cluster_marker_21, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/21_markers.txt", sep="\t")
write.table(cluster_marker_22, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/22_markers.txt", sep="\t")
write.table(cluster_marker_23, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/23_markers.txt", sep="\t")
write.table(cluster_marker_24, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/24_markers.txt", sep="\t")
write.table(cluster_marker_26, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/26_markers.txt", sep="\t")
write.table(cluster_marker_27, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/27_markers.txt", sep="\t")
write.table(cluster_marker_28, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/28_markers.txt", sep="\t")
write.table(cluster_marker_29, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/cluster_markers/29_markers.txt", sep="\t")






#extract the embedding files and metadatafiles for use in other applications
final_metadata <- EA@meta.data
write.table(final_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/Filtered/EA_final_metadata.txt", sep = '\t')

embeddings <- EA@reductions$umap@cell.embeddings
write.table(embeddings, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/Filtered/cell_embedings.txt", sep = '\t')


#umap plot with the adjusted colors. 
DimPlot(object = EA, reduction = "umap", group.by = 'res', cols = c("6" = "#E57373", "10" = "#4CAF50", "26"= "#64B5F6", "27" = "#64B5F6", "28" = "#64B5F6",  "3" = "#FF6F00", "0" = "#3F51B5", "24"= "#BA68C8", "4" = "#4DD0E1","5" = "#4DD0E1", "29" = "#26A69A", "22" = "#F9A825", "16" = "#CDDC39", "8" = "#CDDC39", "20" = "#9E9E9E", "23" = "#9E9E9E", "2" = "#9E9E9E", "19" = "#9E9E9E", "9" = "#9E9E9E", "11" = "#9E9E9E", "1" = "#9E9E9E", "15" = "#9E9E9E", "7" = "#9E9E9E", "13" = "#9E9E9E", "18" = "#9E9E9E", "21" = "#9E9E9E", "17" = "#9E9E9E", "12" = "#9E9E9E"))






###before creating the heatmaps, we need to first make the active.idents for the juvenile dataset to be the column called "res". To do this, we use the SetIdent function. and as a sanity check we see what the active identities are based on juv.data@active.ident. 
EA <- SetIdent(EA, value = EA@meta.data$res)
EA@active.ident


#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(EA, idents ="29")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/Filtered/pearsons_input/EA_29.tmp.txt", sep = '\t')




#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(EA, idents ="29")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/Filtered/pearsons_input/EA_29.tmp.txt", sep = '\t')





###########################################################
#subsetting out neoblasts
EA_neoblast<- subset(EA, idents ="10")

#sanity check neoblast umap
DimPlot(object = EA_neoblast, reduction = "umap", label = T)


#recluster the neoblasts

EA_neoblast <- NormalizeData(object = EA_neoblast, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
EA_neoblast <- FindVariableFeatures(object = EA_neoblast, do.plot = F)

# Scaling the data and removing unwanted sources of variation
EA_neoblast <- ScaleData(object = EA_neoblast, display.progress = F)

# Perform linear dimensional reduction
EA_neoblast <- RunPCA(object = EA_neoblast, pc.genes = juv.data@var.genes, do.print = TRUE)

EA_neoblast <- JackStraw(EA_neoblast, num.replicate = 100)
EA_neoblast <- ScoreJackStraw(EA_neoblast, dims = 1:20)

JackStrawPlot(EA_neoblast, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

EA_neoblast <- FindNeighbors(object = EA_neoblast, dims = c(1:6,8,9,11,13,19))
EA_neoblast <- FindClusters(EA_neoblast, resolution = 0.4, print.output = 0, save.SNN = T)
EA_neoblast <-  RunUMAP(EA_neoblast, reduction.use = "pca", dims = c(1:6,8,9,11,13,19), min.dist = 0.2, n.neighbors = 50)

DimPlot(object = EA_neoblast, reduction = "umap", label = T)

VlnPlot(EA_neoblast, features = "98043523-piwl1-2")
RidgePlot(object = EA_neoblast, features = "98043523-piwl1-2")


#test the markers and their expression
clustermarkers <- FindMarkers(object = EA_neoblast, ident.1 = "1" ,min.pct = 0.25)
FeaturePlot(object=EA_neoblast, c("98050114-h33-3"), pt.size = 0.5)


clusters <- (0:3)
for (i in 0:3) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = EA_neoblast, ident.1 =i, min.pct = 0.25)
  )}

FeaturePlot(object=EA_neoblast, c("98051971-fabp5"), pt.size = 0.5)


write.table(cluster_marker_3, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/neoblast_cluster_markers/3_markers.txt", sep="\t")




#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(EA_neoblast, idents ="0")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/EA/filtered/neoblast_pearsons/EA_0.tmp.txt", sep = '\t')


clustermarkers <- FindMarkers(object = EA_neoblast, ident.1 = "0" ,min.pct = 0.25)
write.table(clustermarkers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/LJ/filtered/neoblast_pearsons/clust_0_markers.txt", sep = '\t')


