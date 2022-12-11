library(Seurat)

x72hpa_matrix <- read.table("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/sc.counts_72hpa.matrix", sep = "\t")


reg_72hpa <- CreateSeuratObject(counts = x72hpa_matrix, min.cells = 5, min.features =200)


# Filter the data
####script for embedding the post embryonic datasets myself
VlnPlot(reg_72hpa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
reg_72hpa <- subset(reg_72hpa, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000)


# Normalizing the data
reg_72hpa <- NormalizeData(object = reg_72hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_72hpa <- FindVariableFeatures(object =reg_72hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_72hpa <- ScaleData(object = reg_72hpa, display.progress = F)

# Perform linear dimensional reduction
reg_72hpa <- RunPCA(object = reg_72hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

reg_72hpa <- JackStraw(reg_72hpa, num.replicate = 100)
reg_72hpa <- ScoreJackStraw(reg_72hpa, dims = 1:20)

JackStrawPlot(reg_72hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

reg_72hpa <- FindNeighbors(object = reg_72hpa, dims = c(1:20))
reg_72hpa <- FindClusters(reg_72hpa, resolution = 1.8, print.output = 0, save.SNN = T)
reg_72hpa <-  RunUMAP(reg_72hpa, reduction.use = "pca", dims= c(1:20), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = reg_72hpa, reduction = "umap", cols = c("26" = "#E57373", "8"= "#E57373" ,"2"= "#E57373" , "3" = "#4CAF50", "21"= "#64B5F6", "22" = "#64B5F6", "5" = "#64B5F6", "16" = "#64B5F6", "11" = "#FF6F00", "7" = "#3F51B5", "25"= "#BA68C8", "0" = "#4DD0E1","14" = "#4DD0E1", "6" = "#26A69A", "10" = "#F9A825", "4" = "#CDDC39","9" = "#CDDC39", "15" = "#CDDC39", "18" = "#9E9E9E", "19" = "#9E9E9E", "17" = "#9E9E9E", "24" = "#9E9E9E", "12" = "#9E9E9E", "1" = "#9E9E9E", "23" = "#9E9E9E", "20" = "#9E9E9E", "27" = "#9E9E9E", "28" = "#9E9E9E", "13" = "#9E9E9E"))


cluster0.markers <- FindMarkers(object = reg_72hpa, ident.1 =0, min.pct = 0.25)

FeaturePlot(reg_72hpa, c("98048987-catl2-4"), pt.size = 0.5)



VlnPlot(reg_72hpa, features = "98043523-piwl1-2")
RidgePlot(object = reg_72hpa, features = "98043523-piwl1-2")




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
neoblast_72hpa<- subset(reg_72hpa, idents=3)
DimPlot(object = neoblast_72hpa, reduction = "umap",label = T)

neoblast_72hpa <- NormalizeData(object = neoblast_72hpa, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
neoblast_72hpa <- FindVariableFeatures(object = neoblast_72hpa, do.plot = F)

# Scaling the data and removing unwanted sources of variation
neoblast_72hpa <- ScaleData(object = neoblast_72hpa, display.progress = F)

# Perform linear dimensional reduction
neoblast_72hpa <- RunPCA(object = neoblast_72hpa, pc.genes = juv.data@var.genes, do.print = TRUE)

neoblast_72hpa <- JackStraw(neoblast_72hpa, num.replicate = 100)
neoblast_72hpa <- ScoreJackStraw(neoblast_72hpa, dims = 1:20)

JackStrawPlot(neoblast_72hpa, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

neoblast_72hpa <- FindNeighbors(object = neoblast_72hpa, dims = c(1:11,13,15,16,20))
neoblast_72hpa <- FindClusters((neoblast_72hpa), resolution = 1.31, print.output = 0, save.SNN = T)
neoblast_72hpa <-  RunUMAP(neoblast_72hpa, reduction.use = "pca", dims=  c(1:11,13,15,16,20), min.dist = 0.3, n.neighbors = 50)

DimPlot(object = neoblast_72hpa, reduction = "umap",label = T)

VlnPlot(neoblast_72hpa, features = "98043523-piwl1-2")
RidgePlot(object = neoblast_72hpa, features = "98043523-piwl1-2")


#testing cluster markers
cluster0.markers <- FindMarkers(object = neoblast_72hpa, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = neoblast_72hpa, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = neoblast_72hpa, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = neoblast_72hpa, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = neoblast_72hpa, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = neoblast_72hpa, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object =neoblast_72hpa, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object =neoblast_72hpa, ident.1 =7, min.pct = 0.25)
cluster8.markers <- FindMarkers(object =neoblast_72hpa, ident.1 =8, min.pct = 0.25)

FeaturePlot(neoblast_72hpa, c("98051971-fabp5"), pt.size = 0.5)


#writing the cluster markers into files
cluster0.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "0", min.pct = 0.25)
write.table(cluster0.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust0_markers.txt", sep="\t")

cluster1.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "1", min.pct = 0.25)
write.table(cluster1.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust1_markers.txt", sep="\t")


cluster2.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "2", min.pct = 0.25)
write.table(cluster2.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust2_markers.txt", sep="\t")


cluster3.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "3", min.pct = 0.25)
write.table(cluster3.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust3_markers.txt", sep="\t")


cluster4.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "4", min.pct = 0.25)
write.table(cluster4.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust4_markers.txt", sep="\t")


cluster5.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "5", min.pct = 0.25)
write.table(cluster5.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust5_markers.txt", sep="\t")


cluster6.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "6", min.pct = 0.25)
write.table(cluster6.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust6_markers.txt", sep="\t")


cluster7.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "7", min.pct = 0.25)
write.table(cluster7.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust7_markers.txt", sep="\t")

cluster8.markers <- FindMarkers(object = neoblast_72hpa, ident.1 = "8", min.pct = 0.25)
write.table(cluster8.markers, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/72hpa/cluster2_subset/cluster_markers/clust8_markers.txt", sep="\t")

#create pearsons correlation metacell files
# Subset the cluster into its own seurat object
cluster_x<- subset(neoblast_72hpa, idents ="7")
# Export the matrix of selected cluster
matrix_cx <- as.matrix(x= cluster_x@assays$RNA@data)
#make sure you keep the file name formatted in the same way as the example below. it should be "abbreviationforyourobject_clusternumber.tmp.txt". This is important for the scripts you use in the next part of the analysis.
#It is also important that you place these files into the same directory as the other scripts provided.
write.table(matrix_cx, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_reg_neoblast/pearsons/x72hpa_7N.tmp.txt", sep = '\t')


