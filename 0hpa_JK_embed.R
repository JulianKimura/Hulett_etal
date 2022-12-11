library(Seurat)

x0hpa_matrix <- read.table("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/0hpa/sc.counts_0hpa.matrix", sep = "\t")


reg_0hpa <- CreateSeuratObject(counts = x0hpa_matrix, min.cells = 5, min.features =200)

VlnPlot(reg_0hpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
reg_0hpa <- subset(reg_0hpa, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 50000)
VlnPlot(reg_0hpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)



# Filter the data
####script for embedding the post embryonic datasets myself


reg_0hpa <- subset(reg_0hpa, subset.names = "nGene", low.thresholds = 500)

# Normalizing the data
reg_0hpa <- NormalizeData(object = reg_0hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_0hpa <- FindVariableFeatures(object = reg_0hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_0hpa <- ScaleData(object = reg_0hpa, display.progress = F)

# Perform linear dimensional reduction
reg_0hpa <- RunPCA(object = reg_0hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

reg_0hpa <- JackStraw(reg_0hpa, num.replicate = 100)
reg_0hpa <- ScoreJackStraw(reg_0hpa, dims = 1:20)

JackStrawPlot(reg_0hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

reg_0hpa <- FindNeighbors(object = reg_0hpa, dims = c(1:20))
reg_0hpa <- FindClusters(reg_0hpa, resolution = 1.5, print.output = 0, save.SNN = T)
reg_0hpa <-  RunUMAP(reg_0hpa, reduction.use = "pca", dims= c(1:20), min.dist = 0.5, n.neighbors = 50)

DimPlot(object = reg_0hpa, reduction = "umap",label = T)


cluster26.markers <- FindMarkers(object = reg_0hpa, ident.1 =26, min.pct = 0.25)

FeaturePlot(reg_0hpa, c("98048987-catl2-4"), pt.size = 0.5)

rm(list=ls(pattern="cluster"))


VlnPlot(reg_0hpa, features = "98043523-piwl1-2")
RidgePlot(object = reg_0hpa, features = "98043523-piwl1-2")


####for loop for making cluster markers and writing them into txt files
clusters <- (0:32)
for (i in 0:32) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = reg_0hpa, ident.1 =i, min.pct = 0.25)
)}

cluster_marker_list <- list()
for (i in 0:32) {append(cluster_marker_list, paste("cluster_marker", i, sep = "_") )}


write.table( cluster_marker_32, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/cluster_markers/clust32_markers.txt", sep="\t")



#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- SubsetData(reg_0hpa, ident.use ="3")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/URD/REH_paper/RH_For_JK_URD_Data2/JK_embedding/neoblast_0hpa/subset_cluster9/neoblast_0hpa_3.tmp.txt", sep = '\t')




########################## SUBSETTING BASED ON NEOBLAST CLUSTER (CLUSTER 0)
#first, subset based on cluster
neoblast_0hpa<- subset(reg_0hpa, idents ="2")
DimPlot(object = neoblast_0hpa, reduction = "umap",label = T)

neoblast_0hpa <- NormalizeData(object = neoblast_0hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
neoblast_0hpa <- FindVariableFeatures(object = neoblast_0hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
neoblast_0hpa <- ScaleData(object = neoblast_0hpa, display.progress = F)

# Perform linear dimensional reduction
neoblast_0hpa <- RunPCA(object = neoblast_0hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

neoblast_0hpa <- JackStraw(neoblast_0hpa, num.replicate = 100)
neoblast_0hpa <- ScoreJackStraw(neoblast_0hpa, dims = 1:20)

JackStrawPlot(neoblast_0hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

neoblast_0hpa <- FindNeighbors(object = neoblast_0hpa, dims = c(1:7,9:11,13,15,18,19))
neoblast_0hpa <- FindClusters(neoblast_0hpa, resolution = 0.6, print.output = 0, save.SNN = T)
neoblast_0hpa <-  RunUMAP(neoblast_0hpa, reduction.use = "pca", dims= c(1:7,9:11,13,15,18,19), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = reg_0hpa, reduction = "umap", cols = c("1" = "#E57373", "9"= "#E57373" ,"32"= "#E57373" ,  "2" = "#4CAF50", "19"= "#64B5F6", "30" = "#64B5F6", "24" = "#64B5F6", "13" = "#64B5F6", "20" = "#FF6F00","4" = "#FF6F00","23" = "#FF6F00", "16" = "#3F51B5", "27"= "#BA68C8", "11" = "#4DD0E1","10" = "#4DD0E1","0" = "#4DD0E1", "22" = "#26A69A", "15" = "#F9A825", "3" = "#CDDC39","8" = "#CDDC39", "31" = "#9E9E9E","22" = "#9E9E9E", "32" = "#9E9E9E", "21" = "#9E9E9E", "18" = "#9E9E9E", "29" = "#9E9E9E", "17" = "#9E9E9E", "28" = "#9E9E9E", "25" = "#9E9E9E", "12" = "#9E9E9E", "5" = "#9E9E9E", "14" = "#9E9E9E", "6" = "#9E9E9E", "7" = "#9E9E9E"))


VlnPlot(neoblast_0hpa, features = "98043523-piwl1-2")
RidgePlot(object = neoblast_0hpa, features = "98043523-piwl1-2")



DimPlot(neoblast_0hpa, reduction = "umap", label = T)
#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object = neoblast_0hpa, ident.1 =7, min.pct = 0.25)


FeaturePlot(neoblast_0hpa, c("98046269-ryr1"), pt.size = 0.5)




#writing the cluster markers into files
cluster0.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "0", min.pct = 0.25)
write.table(cluster0.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/filtered/neoblast/clust0_markers.txt", sep="\t")

cluster1.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "1", min.pct = 0.25)
write.table(cluster1.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/filtered/neoblast/clust1_markers.txt", sep="\t")


cluster2.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "2", min.pct = 0.25)
write.table(cluster2.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/filtered/neoblast/clust2_markers.txt", sep="\t")


cluster3.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "3", min.pct = 0.25)
write.table(cluster3.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/filtered/neoblast/clust3_markers.txt", sep="\t")


cluster4.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "4", min.pct = 0.25)
write.table(cluster4.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/filtered/neoblast/clust4_markers.txt", sep="\t")


cluster5.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "5", min.pct = 0.25)
write.table(cluster5.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/filtered/neoblast/clust5_markers.txt", sep="\t")


cluster6.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "6", min.pct = 0.25)
write.table(cluster6.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/cluster2_susbet/cluster_markers/clust6_markers.txt", sep="\t")


cluster7.markers <- FindMarkers(object = neoblast_0hpa, ident.1 = "7", min.pct = 0.25)
write.table(cluster7.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/cluster2_susbet/cluster_markers/clust7_markers.txt", sep="\t")






#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(neoblast_0hpa, idents ="5")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx,"/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/0hpa/filtered/neoblast/x0hpa_5N.tmp.txt", sep = '\t')


