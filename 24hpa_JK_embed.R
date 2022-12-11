library(Seurat)

x24hpa_matrix <- read.table("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/sc.counts_24hpa.matrix", sep = "\t")


reg_24hpa <- CreateSeuratObject(counts = x24hpa_matrix, min.cells = 5, min.features =200)


# Filter the data
####script for embedding the post embryonic datasets myself
VlnPlot(reg_24hpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
reg_24hpa <- subset(reg_24hpa, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 50000)


# Normalizing the data
reg_24hpa <- NormalizeData(object = reg_24hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_24hpa <- FindVariableFeatures(object = reg_24hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_24hpa <- ScaleData(object = reg_24hpa, display.progress = F)

# Perform linear dimensional reduction
reg_24hpa <- RunPCA(object = reg_24hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

reg_24hpa <- JackStraw(reg_24hpa, num.replicate = 100)
reg_24hpa <- ScoreJackStraw(reg_24hpa, dims = 1:20)

JackStrawPlot(reg_24hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

reg_24hpa <- FindNeighbors(object = reg_24hpa, dims = c(1:20))
reg_24hpa <- FindClusters(reg_24hpa, resolution = 1.2, print.output = 0, save.SNN = T)
reg_24hpa <-  RunUMAP(reg_24hpa, reduction.use = "pca", dims= c(1:20), min.dist = 0.5, n.neighbors = 50)

DimPlot(object = reg_24hpa, reduction = "umap",label = T)


cluster11.markers <- FindMarkers(object = reg_24hpa, ident.1 =11, min.pct = 0.25)

FeaturePlot(reg_24hpa, c("98048987-catl2-4"), pt.size = 0.5)



VlnPlot(reg_24hpa, features = "98043523-piwl1-2")
RidgePlot(object = reg_24hpa, features = "98043523-piwl1-2")




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
neoblast_24hpa<- subset(reg_24hpa, idents =2)
DimPlot(object = neoblast_24hpa, reduction = "umap",label = T)

neoblast_24hpa <- NormalizeData(object = neoblast_24hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
neoblast_24hpa <- FindVariableFeatures(object = neoblast_24hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
neoblast_24hpa <- ScaleData(object = neoblast_24hpa, display.progress = F)

# Perform linear dimensional reduction
neoblast_24hpa <- RunPCA(object = neoblast_24hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

neoblast_24hpa <- JackStraw(neoblast_24hpa, num.replicate = 100)
neoblast_24hpa <- ScoreJackStraw(neoblast_24hpa, dims = 1:20)

JackStrawPlot(neoblast_24hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

neoblast_24hpa <- FindNeighbors(object = neoblast_24hpa, dims = c(1:9,12,14))
neoblast_24hpa <- FindClusters(neoblast_24hpa, resolution = 1.0, print.output = 0, save.SNN = T)
neoblast_24hpa <-  RunUMAP(neoblast_24hpa, reduction.use = "pca", dims= c(1:9,12,14), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = reg_24hpa, reduction = "umap", cols = c("21" = "#E57373", "4"= "#E57373" , "2" = "#4CAF50", "16"= "#64B5F6", "17" = "#64B5F6", "22" = "#64B5F6", "12" = "#64B5F6", "15" = "#FF6F00", "8" = "#FF6F00", "1" = "#3F51B5", "5" = "#4DD0E1", "6"=  "#26A69A",  "14" = "#F9A825", "0" = "#CDDC39", "11" = "#9E9E9E", "9" = "#9E9E9E", "7" = "#9E9E9E", "3" = "#9E9E9E", "24" = "#9E9E9E", "25" = "#9E9E9E", "20" = "#9E9E9E", "23" = "#9E9E9E", "10" = "#9E9E9E", "13" = "#9E9E9E", "18" = "#9E9E9E", "19" = "#9E9E9E"))

VlnPlot(neoblast_24hpa, features = "98043523-piwl1-2")
RidgePlot(object = neoblast_24hpa, features = "98043523-piwl1-2")


#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =7, min.pct = 0.25)
cluster8.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =8, min.pct = 0.25)

FeaturePlot(neoblast_24hpa, c("98005177-actb"), pt.size = 0.5)


#writing the cluster markers into files
#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object = neoblast_24hpa, ident.1 =7, min.pct = 0.25)

FeaturePlot(neoblast_24hpa, c("98050114-h33-3"), pt.size = 0.5)




#writing the cluster markers into files
cluster0.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "0", min.pct = 0.25)
write.table(cluster0.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust0_markers.txt", sep="\t")

cluster1.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "1", min.pct = 0.25)
write.table(cluster1.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust1_markers.txt", sep="\t")


cluster2.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "2", min.pct = 0.25)
write.table(cluster2.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust2_markers.txt", sep="\t")


cluster3.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "3", min.pct = 0.25)
write.table(cluster3.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust3_markers.txt", sep="\t")


cluster4.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "4", min.pct = 0.25)
write.table(cluster4.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust4_markers.txt", sep="\t")


cluster5.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "5", min.pct = 0.25)
write.table(cluster5.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust5_markers.txt", sep="\t")


cluster6.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "6", min.pct = 0.25)
write.table(cluster6.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust6_markers.txt", sep="\t")


cluster7.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "7", min.pct = 0.25)
write.table(cluster7.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust7_markers.txt", sep="\t")

cluster8.markers <- FindMarkers(object = neoblast_24hpa, ident.1 = "8", min.pct = 0.25)
write.table(cluster8.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/24hpa/cluster4_subset/cluster_markers/clust8_markers.txt", sep="\t")



#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(neoblast_24hpa, idents ="5")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_reg_neoblast/pearsons/x24hpa_5N.tmp.txt", sep = '\t')
