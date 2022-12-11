# Load libraries
library(Seurat)
library(magrittr)
library(dplyr)
require(scales)


# Load the dataset
juv_matrix.data <- read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/sc.counts_J.matrix", sep = "\t")

juv.data <- CreateSeuratObject(counts = juv_matrix.data, min.cells = 5, min.features = 200)
VlnPlot(juv.data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)


# Filter the data. the vlnplot shows that there are some outliers where the cells have more than 5000 genes and more than 15000 reads. these suggests doublets, or two cells that were sequenced togehter. We also see that the majority of the cells have more than 500 genes. We discard any cells with less than 500 because these were most likely to be low quality samples, which would add noise to our data. 
juv.data <- subset(juv.data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 15000)


# Normalizing the data
juv.data <- NormalizeData(object = juv.data, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
juv.data <- FindVariableFeatures(object = juv.data, do.plot = F)

# Scaling the data and removing unwanted sources of variation
juv.data <- ScaleData(object = juv.data, display.progress = F)

# Perform linear dimensional reduction
juv.data <- RunPCA(object = juv.data, pc.genes = juv.data@var.genes, do.print = TRUE)


#jackstraw
juv.data <- JackStraw(juv.data, num.replicate = 100)
juv.data <- ScoreJackStraw(juv.data, dims = 1:20)
JackStrawPlot(juv.data, dims = 1:20)

#calc umap
juv.data <- FindNeighbors(object = juv.data, dims = 1:20)
juv.data <- FindClusters(juv.data, resolution = 1.2, print.output = 0, save.SNN = T)
juv.data <-  RunUMAP(juv.data, reduction.use = "pca", dims= 1:20, min.dist = 0.6, n.neighbors = 50)

DimPlot(object = juv.data, reduction = "umap", label = T, group.by = 'res')

#test the markers and their expression
clustermarkers <- FindMarkers(object = juv.data, ident.1 = "21" ,min.pct = 0.25)
FeaturePlot(object=juv.data, c("98046213-foxo1"), pt.size = 0.5)


#I found that at the resolution parameter of 0.5, we are underclustering a lot of the data, but we have the gut cluster being retained as whole. We will be using this dataset with res parameter of 0.5 as the base for making a "custom" metadata file. At resolution parameter of 1.2, we see other clusters being 

#first, extract the metadata containing the clustering information at resolution of 0.5
res0.5_metadata <- juv.data@meta.data
#now I am trimming the metadata so that there is only the clustering information for the resolution of 0.5 present.
res0.5_metadata <- res0.5_metadata[,-5]
#write this metadatafile in the desired directory. this will be concatentated with clustering information from other clustering parameters using a python script. 
write.table(res0.5_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/res0.5_metadata.txt", sep="\t")


#Next, I want to subset the dataset with the resolution parameter of 1.2. The idea is to subset out the clusters that I want to add-in to the 0.5 resolution UMAP. In the 1.2 
subsetted_clusters<- subset(juv.data, idents = c("24", "22"))
#plot the subsetted data on a umap for a sanity check. 
DimPlot(object = subsetted_clusters, reduction = "umap", label = T)

#now we want to extract the metadata for these two subsetted clusters
res1.2_metadata <- subsetted_clusters@meta.data
#now I am trimming the metadata so that there is only the clustering information for the resolution of 0.5 present.
res1.2_metadata <- res1.2_metadata[,-5]

#Next, I want to rename the cluster identities here to retain a level of continuity. at a resolution parameter of 0.5, we have 19 clusters. So we want to rename cluster 24 as 20, and cluster 22 as 21. 
current.cluster.ids <- c("24", "22")
new.cluster.ids <- c("20", "21")
res1.2_metadata$RNA_snn_res.1.2<- plyr::mapvalues(x = res1.2_metadata$RNA_snn_res.1.2, from = current.cluster.ids, to = new.cluster.ids)


#write this metadatafile in the desired directory. this will be concatentated with clustering information from other clustering parameters using a python script. 
write.table(res1.2_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/res1.2_metadata.txt", sep="\t")


#now that the custom metadata file has been generated, we re-import this back into R. 
custom_metadata <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/custom_metadata.txt", sep = "\t")


#extract the column from the custom_metadata file with the clustering information, and rename it res
res <- custom_metadata[,4]

unique(res)
#place this column into the juv.data seurat object's metadata so you can call it when plotting a UMAP. 
juv.data@meta.data <- cbind(juv.data@meta.data, res)


#now generate the custom UMAP plot with the custom clustering information. 
DimPlot(object = juv.data, reduction = "umap", label = T, group.by = 'res', cols = c("0" = "#F8766D", "1" =  "#EB8335", "2" = "#DA8F00", "3" = "#C49A00", "4" = "#A9A400", "5" = "#86AC00", "6" = "#53B400", "7" = "#00BA38", "8" = "#00BE6D", "9" = "#00C094", "10" =  "#00C0B5", "11"= "#00BDD2", "12" =  "#00B6EB", "13" =  "#00ABFD", "14" = "#619CFF", "18" = "#9E9E9E", "19" = "#9E9E9E", "20" = "#9E9E9E", "21" = "#9E9E9E", "15" = "#9E9E9E", "17" = "#9E9E9E", "16" = "#9E9E9E"))

library(scales)
show_col(hue_pal()(21))


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=21)
color_list




#extract the embedding files and metadatafiles for use in other applications
final_metadata <- juv.data@meta.data
write.table(final_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/HJ_final_metadata.txt", sep = '\t')

embeddings <- juv.data@reductions$umap@cell.embeddings
write.table(embeddings, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cell_embedings.txt", sep = '\t')

hj_count_matrix <- juv.data@assays$RNA@counts
write.table(hj_count_matrix, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/hj_count_matrix.txt", sep = '\t')


#umap plot with the adjusted colors. 
DimPlot(object = juv.data, reduction = "umap", group.by = 'res', cols = c("1" = "#E57373", "12"= "#E57373" , "3" = "#4CAF50", "7"= "#64B5F6", "8" = "#64B5F6", "6" = "#64B5F6", "2" = "#64B5F6", "11" = "#FF6F00", "5" = "#3F51B5", "13"= "#BA68C8", "14" = "#4DD0E1", "10" = "#26A69A", "9" = "#BA68C8", "4" = "#F9A825", "0" = "#CDDC39", "18" = "#9E9E9E", "19" = "#9E9E9E", "20" = "#9E9E9E", "21" = "#9E9E9E", "15" = "#9E9E9E", "17" = "#9E9E9E", "16" = "#9E9E9E"))




###before creating the heatmaps, we need to first make the active.idents for the juvenile dataset to be the column called "res". To do this, we use the SetIdent function. and as a sanity check we see what the active identities are based on juv.data@active.ident. 
juv.data <- SetIdent(juv.data, value = juv.data@meta.data$res)
juv.data@active.ident


#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(juv.data, idents ="21")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/pearsons_input/HJ_21.tmp.txt", sep = '\t')


#for loop for finding markers
clusters <- (0:21)
for (i in 0:21) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = juv.data, ident.1 =i, min.pct = 0.25)
  )}



write.table(cluster_marker_0, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/0_markers.txt", sep="\t")
write.table(cluster_marker_1, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/1_markers.txt", sep="\t")
write.table(cluster_marker_2, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/2_markers.txt", sep="\t")
write.table(cluster_marker_3, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/3_markers.txt", sep="\t")
write.table(cluster_marker_4, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/4_markers.txt", sep="\t")
write.table(cluster_marker_5, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/5_markers.txt", sep="\t")
write.table(cluster_marker_6, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/6_markers.txt", sep="\t")
write.table(cluster_marker_7, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/7_markers.txt", sep="\t")
write.table(cluster_marker_8, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/8_markers.txt", sep="\t")
write.table(cluster_marker_9, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/9_markers.txt", sep="\t")
write.table(cluster_marker_10, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/10_markers.txt", sep="\t")
write.table(cluster_marker_11, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/11_markers.txt", sep="\t")
write.table(cluster_marker_12, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/12_markers.txt", sep="\t")
write.table(cluster_marker_13, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/13_markers.txt", sep="\t")
write.table(cluster_marker_14, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/14_markers.txt", sep="\t")
write.table(cluster_marker_15, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/15_markers.txt", sep="\t")
write.table(cluster_marker_16, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/16_markers.txt", sep="\t")
write.table(cluster_marker_17, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/17_markers.txt", sep="\t")
write.table(cluster_marker_18, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/18_markers.txt", sep="\t")
write.table(cluster_marker_19, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/19_markers.txt", sep="\t")
write.table(cluster_marker_20, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/20_markers.txt", sep="\t")
write.table(cluster_marker_21, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/cluster_markers/21_markers.txt", sep="\t")





###########################################################
#subsetting out neoblasts
HJ_neoblast<- subset(juv.data, idents ="3")

#sanity check neoblast umap
DimPlot(object = HJ_neoblast, reduction = "umap", label = T)


#recluster the neoblasts

HJ_neoblast <- NormalizeData(object = HJ_neoblast, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
HJ_neoblast <- FindVariableFeatures(object = HJ_neoblast, do.plot = F)

# Scaling the data and removing unwanted sources of variation
HJ_neoblast <- ScaleData(object = HJ_neoblast, display.progress = F)

# Perform linear dimensional reduction
HJ_neoblast <- RunPCA(object = HJ_neoblast, pc.genes = juv.data@var.genes, do.print = TRUE)

HJ_neoblast <- JackStraw(HJ_neoblast, num.replicate = 100)
HJ_neoblast <- ScoreJackStraw(HJ_neoblast, dims = 1:20)

JackStrawPlot(HJ_neoblast, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

HJ_neoblast <- FindNeighbors(object = HJ_neoblast, dims = c(1:9,11:14,17,18))
HJ_neoblast <- FindClusters(HJ_neoblast, resolution = 0.9, print.output = 0, save.SNN = T)
HJ_neoblast <-  RunUMAP(HJ_neoblast, reduction.use = "pca", dims= c(1:9,11:14,17,18), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = HJ_neoblast, reduction = "umap", label = T)

VlnPlot(HJ_neoblast, features = "98050813-ccnb2")
RidgePlot(object = HJ_neoblast, features = "98043523-piwl1-2", group.by = 'res')


#test the markers and their expression
FeaturePlot(object=int_neoblast, c("98036104-boll"), pt.size = 0.5)







###############################
##we need to do another round of combining clusters into custom clusters
#here i am setting the baseline as res0.5. I will be adding in clustering from res0.9 into this umap
HJ_neoblast@meta.data <- HJ_neoblast@meta.data[-4]
neoblast_res0.5_metadata <- HJ_neoblast@meta.data
write.table(neoblast_res0.5_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/Neoblast_subset/neoblast_res_0.5_metadata.txt", sep="\t")


#make sure that the res parameter that we are extracting from is 0.9 from this point forward. 
HJ_neoblast <- SetIdent(HJ_neoblast, value = HJ_neoblast@meta.data$RNA_snn_res.0.9)
HJ_neoblast@active.ident


#now I need to subset the clusters from res parameter 0.9. specifically, clusters 1,3,2,5. 
subsetted_clusters<- subset(HJ_neoblast, idents = c("1", "2", "3", "5"))
#plot the subsetted data on a umap for a sanity check. 
DimPlot(object = subsetted_clusters, reduction = "umap", label = T)

#now we want to extract the metadata for these two subsetted clusters
res0.9_metadata <- subsetted_clusters@meta.data
#now I am trimming the metadata so that there is only the clustering information for the resolution of 0.5 present.
res0.9_metadata <- res0.9_metadata[,-5]

#Next, I want to rename the cluster identities here to retain a level of continuity. at a resolution parameter of 0.5, we have 19 clusters. So we want to rename cluster 24 as 20, and cluster 22 as 21. 
current.cluster.ids <- c("1", "3", "2", "5")
new.cluster.ids <- c("1", "2", "4", "5" )
res0.9_metadata$RNA_snn_res.0.9<- plyr::mapvalues(x = res0.9_metadata$RNA_snn_res.0.9, from = current.cluster.ids, to = new.cluster.ids)


#write this metadatafile in the desired directory. this will be concatentated with clustering information from other clustering parameters using a python script. 
write.table(res0.9_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/Neoblast_subset/res0.9_metadata.txt", sep="\t")


#Read in the custom metadata that was generated. 
#now that the custom metadata file has been generated, we re-import this back into R. 
custom_metadata <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/Neoblast_subset/custom_metadata.txt", sep = "\t")


#extract the column from the custom_metadata file with the clustering information, and rename it res
res <- custom_metadata[,4]

unique(res)
#place this column into the juv.data seurat object's metadata so you can call it when plotting a UMAP. 
HJ_neoblast@meta.data <- cbind(HJ_neoblast@meta.data, res)

#######FOR CARLOS#################
#now generate the custom UMAP plot with the custom clustering information. 
DimPlot(object = HJ_neoblast, reduction = "umap", label = T, group.by = 'res')


#make sure the active idents are from res in the neoblast subset. 
HJ_neoblast <- SetIdent(HJ_neoblast, value = HJ_neoblast@meta.data$res)
HJ_neoblast@active.ident

#finding markers
clustermarkers <- FindMarkers(object = HJ_neoblast, ident.1 = "3" ,min.pct = 0.25)
FeaturePlot(object=HJ_neoblast, c("98050813-ccnb2"), pt.size = 0.5)
write.table(clustermarkers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/Neoblast_subset/custom_cluster_markers/clust5.txt", sep = '\t')
###########STOP CARLOS############


clusters <- (0:5)
for (i in 0:5) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = HJ_neoblast, ident.1 =i, min.pct = 0.25)
  )}

write.table(cluster_marker_5, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/Neoblast_subset/pearsons_neoblast_2022/cluster_markers/5_markers.txt", sep = '\t')

write.table(A_neoblast@assays$RNA@counts, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Postembryonic_neoblast_integration_2022/A_neoblast_matrix.txt", sep = '\t')


#to make GO background terms for HJ neoblasts, I am taking the neoblast cluster markers, and then using them as background. 
clustermarkers <- FindMarkers(object = juv.data, ident.1 = "3" ,min.pct = 0.25)
write.table(clustermarkers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/Neoblast_subset/custom_cluster_markers/neoblast_general.txt", sep = '\t')


#creating input files for pearsons correlation
#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(HJ_neoblast, idents ="5")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/Neoblast_subset/pearsons_neoblast_2022/HJ_5.tmp.txt", sep = '\t')

