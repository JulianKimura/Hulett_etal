library(Seurat)

x6hpa_matrix <- read.table("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/6hpa/sc.counts_6hpa.matrix", sep = "\t")


reg_6hpa <- CreateSeuratObject(counts = x6hpa_matrix, min.cells = 5, min.features =200)


# Filter the data
####script for embedding the post embryonic datasets myself
VlnPlot(reg_6hpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
reg_6hpa <- subset(reg_6hpa, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 40000)

# Normalizing the data
reg_6hpa <- NormalizeData(object = reg_6hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_6hpa <- FindVariableFeatures(object = reg_6hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_6hpa <- ScaleData(object = reg_6hpa, display.progress = F)

# Perform linear dimensional reduction
reg_6hpa <- RunPCA(object = reg_6hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

reg_6hpa <- JackStraw(reg_6hpa, num.replicate = 100)
reg_6hpa <- ScoreJackStraw(reg_6hpa, dims = 1:20)

JackStrawPlot(reg_6hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

reg_6hpa <- FindNeighbors(object = reg_6hpa, dims = c(1:20))
reg_6hpa <- FindClusters(reg_6hpa, resolution = 1.7, print.output = 0, save.SNN = T)
reg_6hpa <-  RunUMAP(reg_6hpa, reduction.use = "pca", dims= c(1:20), min.dist = 0.5, n.neighbors = 50)

DimPlot(reg_6hpa, reduction = 'umap')

DimPlot(object = reg_6hpa, reduction = "umap", cols = c("5" = "#E57373", "15"= "#E57373" , "4" = "#4CAF50", "27"= "#64B5F6", "6" = "#64B5F6", "30" = "#64B5F6", "21" = "#64B5F6", "6" = "#FF6F00","3" = "#FF6F00", "0" = "#3F51B5", "28"= "#BA68C8", "9" = "#4DD0E1", "31" = "#26A69A", "13" = "#F9A825", "2" = "#CDDC39", "20" = "#CDDC39", "32" = "#9E9E9E", "29" = "#9E9E9E", "11" = "#9E9E9E", "22" = "#9E9E9E", "19" = "#9E9E9E", "16" = "#9E9E9E", "1" = "#9E9E9E",  "14" = "#9E9E9E",  "26" = "#9E9E9E",  "10" = "#9E9E9E",  "7" = "#9E9E9E",  "18" = "#9E9E9E",  "25" = "#9E9E9E",  "12" = "#9E9E9E",  "17" = "#9E9E9E",  "23" = "#9E9E9E",  "24" = "#9E9E9E"))


rm(list=ls(pattern="cluster"))


cluster16.markers <- FindMarkers(object = reg_6hpa, ident.1 =16, min.pct = 0.25)

FeaturePlot(reg_6hpa, c("98048987-catl2-4"), pt.size = 0.5)



VlnPlot(reg_6hpa, features = "98043523-piwl1-2")
RidgePlot(object = reg_6hpa, features = "98043523-piwl1-2")



####for loop for making cluster markers and writing them into txt files
for (i in 0:32) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = reg_6hpa, ident.1 =i, min.pct = 0.25)
  )}

write.table(cluster_marker_0, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/0_markers.txt", sep="\t")
write.table(cluster_marker_1, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/1_markers.txt", sep="\t")
write.table(cluster_marker_2, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/2_markers.txt", sep="\t")
write.table(cluster_marker_3, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/3_markers.txt", sep="\t")
write.table(cluster_marker_4, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/4_markers.txt", sep="\t")
write.table(cluster_marker_5, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/5_markers.txt", sep="\t")
write.table(cluster_marker_6, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/6_markers.txt", sep="\t")
write.table(cluster_marker_7, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/7_markers.txt", sep="\t")
write.table(cluster_marker_8, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/8_markers.txt", sep="\t")
write.table(cluster_marker_9, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/9_markers.txt", sep="\t")
write.table(cluster_marker_10, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/10_markers.txt", sep="\t")
write.table(cluster_marker_11, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/11_markers.txt", sep="\t")
write.table(cluster_marker_12, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/12_markers.txt", sep="\t")
write.table(cluster_marker_13, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/13_markers.txt", sep="\t")
write.table(cluster_marker_14, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/14_markers.txt", sep="\t")
write.table(cluster_marker_15, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/15_markers.txt", sep="\t")
write.table(cluster_marker_16, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/16_markers.txt", sep="\t")
write.table(cluster_marker_17, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/17_markers.txt", sep="\t")
write.table(cluster_marker_18, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/18_markers.txt", sep="\t")
write.table(cluster_marker_19, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/19_markers.txt", sep="\t")
write.table(cluster_marker_20, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/20_markers.txt", sep="\t")
write.table(cluster_marker_21, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/21_markers.txt", sep="\t")
write.table(cluster_marker_22, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/22_markers.txt", sep="\t")
write.table(cluster_marker_23, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/23_markers.txt", sep="\t")
write.table(cluster_marker_24, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/24_markers.txt", sep="\t")
write.table(cluster_marker_25, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/25_markers.txt", sep="\t")
write.table(cluster_marker_26, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/26_markers.txt", sep="\t")
write.table(cluster_marker_27, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/27_markers.txt", sep="\t")
write.table(cluster_marker_28, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/28_markers.txt", sep="\t")
write.table(cluster_marker_29, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/29_markers.txt", sep="\t")
write.table(cluster_marker_30, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/30_markers.txt", sep="\t")
write.table(cluster_marker_31, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/31_markers.txt", sep="\t")
write.table(cluster_marker_32, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/filtered/cluster_markers/32_markers.txt", sep="\t")


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
neoblast_6hpa<- subset(reg_6hpa, idents ="4")
DimPlot(object = neoblast_6hpa, reduction = "umap",label = T)

neoblast_6hpa <- NormalizeData(object = neoblast_6hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
neoblast_6hpa <- FindVariableFeatures(object = neoblast_6hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
neoblast_6hpa <- ScaleData(object = neoblast_6hpa, display.progress = F)

# Perform linear dimensional reduction
neoblast_6hpa <- RunPCA(object = neoblast_6hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

neoblast_6hpa <- JackStraw(neoblast_6hpa, num.replicate = 100)
neoblast_6hpa <- ScoreJackStraw(neoblast_6hpa, dims = 1:20)

JackStrawPlot(neoblast_6hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

neoblast_6hpa <- FindNeighbors(object = neoblast_6hpa, dims = c(1:4,6:8,11,12))
neoblast_6hpa <- FindClusters(neoblast_6hpa, resolution = 0.8, print.output = 0, save.SNN = T)
neoblast_6hpa <-  RunUMAP(neoblast_6hpa, reduction.use = "pca", dims= c(1:4,6:8,11,12), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = neoblast_6hpa, reduction = "umap",label = T)

VlnPlot(neoblast_6hpa, features = "98043523-piwl1-2")
RidgePlot(object = neoblast_6hpa, features = "98043523-piwl1-2")


#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object = neoblast_6hpa, ident.1 =7, min.pct = 0.25)

FeaturePlot(neoblast_6hpa, c("98049625-nfic"), pt.size = 0.5)


#writing the cluster markers into files
cluster0.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "0", min.pct = 0.25)
write.table(cluster0.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust0_markers.txt", sep="\t")

cluster1.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "1", min.pct = 0.25)
write.table(cluster1.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust1_markers.txt", sep="\t")


cluster2.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "2", min.pct = 0.25)
write.table(cluster2.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust2_markers.txt", sep="\t")


cluster3.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "3", min.pct = 0.25)
write.table(cluster3.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust3_markers.txt", sep="\t")


cluster4.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "4", min.pct = 0.25)
write.table(cluster4.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust4_markers.txt", sep="\t")


cluster5.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "5", min.pct = 0.25)
write.table(cluster5.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust5_markers.txt", sep="\t")


cluster6.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "6", min.pct = 0.25)
write.table(cluster6.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust6_markers.txt", sep="\t")


cluster7.markers <- FindMarkers(object = neoblast_6hpa, ident.1 = "7", min.pct = 0.25)
write.table(cluster7.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/6hpa/cluster1_subset/cluster_markers/clust7_markers.txt", sep="\t")




#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(neoblast_6hpa, idents ="4")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_reg_neoblast/pearsons/x6hpa_4N.tmp.txt", sep = '\t')

