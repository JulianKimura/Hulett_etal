x0hpa_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/0hpa_neoblast_matrix", sep = "\t"))
x24hpa_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/24hpa_neoblast_matrix", sep = "\t"))
x6hpa_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/6hpa_neoblast_matrix", sep = "\t"))
x8dpa_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/8dpa_neoblast_matrix", sep = "\t"))
x17dpa_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/17dpa_neoblast_matrix", sep = "\t"))
x29dpa_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/29dpa_neoblast_matrix", sep = "\t"))
x72dpa_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/72hpa_neoblast_matrix", sep = "\t"))

#Merge the datasets
list_stages <- list(x0hpa_matrix,x24hpa_matrix,x6hpa_matrix,x8dpa_matrix,x17dpa_matrix,x29dpa_matrix,x72dpa_matrix)
int_data_matrix <- dsCombineDGE(list_stages)
write.table(int_data_matrix, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/int_data_matrix.txt", sep="\t")

int_data_matrix <- as.matrix(read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Regeneration/filtered_integrated_neoblast/int_data_matrix.txt", sep = "\t"))

reg_filtered.merged <- CreateSeuratObject(counts = int_data_matrix, min.cells = 5, min.features =200)


#reg_filtered.merged<- subset(reg_filtered.merged, subset = nFeature_RNA > 500 & nFeature_RNA < 4000)


# Normalizing the data
reg_filtered.merged <- NormalizeData(object = reg_filtered.merged, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
reg_filtered.merged <- FindVariableFeatures(object = reg_filtered.merged, do.plot = F)

# Scaling the data and removing unwanted sources of variation
reg_filtered.merged <- ScaleData(object = reg_filtered.merged, display.progress = F)

# Perform linear dimensional reduction
reg_filtered.merged <- RunPCA(object = reg_filtered.merged, pc.genes = int_embryo@var.genes, do.print = TRUE)


#jackstraw
reg_filtered.merged <- JackStraw(object = reg_filtered.merged, num.replicate = 100)
reg_filtered.merged <- ScoreJackStraw(object = reg_filtered.merged, dims = 1:20)
JackStrawPlot(object = reg_filtered.merged, dims = 1:20)


#umap
reg_filtered.merged <- FindNeighbors(object = reg_filtered.merged, dims = 1:20)
reg_filtered.merged <- FindClusters(reg_filtered.merged, resolution = 1.4, print.output = 0, save.SNN = T)
reg_filtered.merged <-  RunUMAP(reg_filtered.merged, reduction.use = "pca", dims= 1:20, min.dist= 0.3, n.neighbors =49 )
DimPlot(object = reg_filtered.merged, reduction = "umap",label = T)


FeaturePlot(object=reg_filtered.merged, c("98048987-catl2-4"), pt.size = 0.5)


clusters <- (0:37)
for (i in 0:37) {
  g <- paste("cluster_marker", i, sep = "_")
  assign(g, FindMarkers(object = reg_filtered.merged, ident.1 =i, min.pct = 0.25)
  )}


write.table(cluster_marker_0, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/0_markers.txt", sep = '\t')
write.table(cluster_marker_1, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/1_markers.txt", sep = '\t')
write.table(cluster_marker_2, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/2_markers.txt", sep = '\t')
write.table(cluster_marker_3, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/3_markers.txt", sep = '\t')
write.table(cluster_marker_4, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/4_markers.txt", sep = '\t')
write.table(cluster_marker_5, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/5_markers.txt", sep = '\t')
write.table(cluster_marker_6, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/6_markers.txt", sep = '\t'),
write.table(cluster_marker_7, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/7_markers.txt", sep = '\t')
write.table(cluster_marker_8, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/8_markers.txt", sep = '\t')
write.table(cluster_marker_9, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/9_markers.txt", sep = '\t')
write.table(cluster_marker_10, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/10_markers.txt", sep = '\t')
write.table(cluster_marker_11, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/11_markers.txt", sep = '\t')
write.table(cluster_marker_12, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/12_markers.txt", sep = '\t')
write.table(cluster_marker_13, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/13_markers.txt", sep = '\t')
write.table(cluster_marker_14, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/14_markers.txt", sep = '\t')
write.table(cluster_marker_15, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/15_markers.txt", sep = '\t')
write.table(cluster_marker_16, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/16_markers.txt", sep = '\t')
write.table(cluster_marker_17, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/17_markers.txt", sep = '\t')
write.table(cluster_marker_18, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/18_markers.txt", sep = '\t')
write.table(cluster_marker_19, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/19_markers.txt", sep = '\t')
write.table(cluster_marker_20, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/20_markers.txt", sep = '\t')
write.table(cluster_marker_21, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/21_markers.txt", sep = '\t')
write.table(cluster_marker_22, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/22_markers.txt", sep = '\t')
write.table(cluster_marker_23, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/23_markers.txt", sep = '\t')
write.table(cluster_marker_24, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/24_markers.txt", sep = '\t')
write.table(cluster_marker_25, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/25_markers.txt", sep = '\t')
write.table(cluster_marker_26, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/26_markers.txt", sep = '\t')
write.table(cluster_marker_27, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/27_markers.txt", sep = '\t')
write.table(cluster_marker_28, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/28_markers.txt", sep = '\t')
write.table(cluster_marker_29, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/29_markers.txt", sep = '\t')
write.table(cluster_marker_30, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/30_markers.txt", sep = '\t')
write.table(cluster_marker_31, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/31_markers.txt", sep = '\t')
write.table(cluster_marker_32, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/32_markers.txt", sep = '\t')
write.table(cluster_marker_33, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/33_markers.txt", sep = '\t')
write.table(cluster_marker_34, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/34_markers.txt", sep = '\t')
write.table(cluster_marker_35, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/35_markers.txt", sep = '\t')
write.table(cluster_marker_36, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/36_markers.txt", sep = '\t')
write.table(cluster_marker_37, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/cluster_markers/37_markers.txt", sep = '\t')






DimPlot(object = reg_filtered.merged, reduction = "umap", group.by= "orig.ident")
DimPlot(object = reg_filtered.merged, reduction = "umap", group.by= "orig.ident", split.by = "orig.ident")
DimPlot(object = reg_filtered.merged, reduction = "umap",label = T)


clustermarkers <- FindMarkers(object = reg_filtered.merged, ident.1 = "4" ,min.pct = 0.25)



RidgePlot(object = reg_filtered.merged, features = "98043523-piwl1-2")
VlnPlot(reg_filtered.merged, features = "98043523-piwl1-2")


FeaturePlot(object=reg_filtered.merged, c("98132150-h10-6"), pt.size = 0.5)



#umap plot with the adjusted colors. 
DimPlot(object = reg_filtered.merged, reduction = "umap", cols = c("1" = "#E57373", "9"= "#E57373" , "2" = "#4CAF50", "12"= "#64B5F6", "11" = "#64B5F6", "25" = "#64B5F6", "17" = "#64B5F6", "14" = "#FF6F00",  "8" = "#FF6F00", "6" = "#3F51B5","33" = "#3F51B5", "21"= "#BA68C8", "0" = "#4DD0E1", "7" = "#4DD0E1", "16" = "#26A69A", "20" = "#F9A825","29" = "#F9A825", "10" = "#CDDC39", "19" = "#CDDC39","3" = "#CDDC39","36" = "#CDDC39","34" = "#9E9E9E", "22" = "#9E9E9E", "23" = "#9E9E9E", "28" = "#9E9E9E", "24" = "#9E9E9E", "27" = "#9E9E9E", "31" = "#9E9E9E", "30" = "#9E9E9E", "37" = "#9E9E9E", "18" = "#9E9E9E", "5" = "#9E9E9E", "4" = "#9E9E9E", "32" = "#9E9E9E", "26" = "#9E9E9E", "13" = "#9E9E9E", "35" = "#9E9E9E", "15" = "#9E9E9E"))


#first, subset based on cluster
int_reg_neoblast<- subset(reg_filtered.merged, idents ="2")
DimPlot(object = int_reg_neoblast, reduction = "umap",label = T)

int_reg_neoblast <- NormalizeData(object = int_reg_neoblast, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
int_reg_neoblast <- FindVariableFeatures(object = int_reg_neoblast, do.plot = F)

# Scaling the data and removing unwanted sources of variation
int_reg_neoblast <- ScaleData(object = int_reg_neoblast, display.progress = F)

# Perform linear dimensional reduction
int_reg_neoblast <- RunPCA(object = int_reg_neoblast, pc.genes = juv.data@var.genes, do.print = TRUE)

int_reg_neoblast <- JackStraw(int_reg_neoblast, num.replicate = 100)
int_reg_neoblast <- ScoreJackStraw(int_reg_neoblast, dims = 1:20)

JackStrawPlot(int_reg_neoblast, dims = 1:20)

# Run non-linear dimensional reduction (UMAP)

int_reg_neoblast <- FindNeighbors(object = int_reg_neoblast, dims = c(1:20))
int_reg_neoblast <- FindClusters(int_reg_neoblast, resolution = 0.6, print.output = 0, save.SNN = T)
int_reg_neoblast <-  RunUMAP(int_reg_neoblast, reduction.use = "pca", dims= c(1:20), min.dist = 0.3, n.neighbors = 50)

current.stage.ids <- c("X0hpa1",  "X0hpa2",  "X0hpa3",  "X6hpa1",  "X6hpa2",  "X6hpa3",  "X24hpa1", "X24hpa2", "X24hpa3", "X72hpa1", "X72hpa2", "X72hpa3", "X8dpa1",  "X8dpa2",  "X8dpa3",  "X17dpa1", "X17dpa2", "X17dpa3", "X29dpa1", "X29dpa2", "X29dpa3")
new.stage.ids <- c("X0hpa",  "X0hpa",  "X0hpa",  "X6hpa",  "X6hpa",  "X6hpa",  "X24hpa", "X24hpa", "X24hpa", "X72hpa", "X72hpa", "X72hpa", "X8dpa",  "X8dpa",  "X8dpa",  "X17dpa", "X17dpa", "X17dpa", "X29dpa", "X29dpa", "X29dpa")
int_reg_neoblast@meta.data$orig.ident<- plyr::mapvalues(x = int_reg_neoblast@meta.data$orig.ident, from = current.stage.ids, to = new.stage.ids)



DimPlot(object = int_reg_neoblast, reduction = "umap",label = T)
DimPlot(object = int_reg_neoblast, reduction = "umap",group.by = "orig.ident")

FeaturePlot(int_reg_neoblast, c("98039600-tba4a"), pt.size = 0.5)
FeaturePlot(int_reg_neoblast, c("98050114-h33-3"), pt.size = 0.5)


cluster0.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =0, min.pct = 0.25)
cluster1.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =1, min.pct = 0.25)
cluster2.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =2, min.pct = 0.25)
cluster3.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =3, min.pct = 0.25)
cluster4.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =4, min.pct = 0.25)
cluster5.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =5, min.pct = 0.25)
cluster6.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =6, min.pct = 0.25)
cluster7.markers <- FindMarkers(object = int_reg_neoblast, ident.1 =7, min.pct = 0.25)

VlnPlot(int_reg_neoblast, features = "98043523-piwl1-2")
RidgePlot(object = int_reg_neoblast, features = "98043523-piwl1-2")


###############I am finding that at higher resolution parameters, we are able to extract a cluster identity associated with boule or the germline. to extract this cluster, we need to first subset the res0.9 parameter first, so it can be added into the res0.6

#First I need to extract the metadata file for the entire dataset when clustered with a resolution of 0.6


res0.6_metadata <- int_reg_neoblast@meta.data
res0.6_metadata <- res0.6_metadata[,-5]

#write this metadata file into the directory of choice
write.table(res0.6_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/res0.6_metadata.txt", sep="\t")

#now subset out the boule cluster identities
boule_subsetted<- subset(int_reg_neoblast, idents = c("10"))
DimPlot(object = boule_subsetted, reduction = "umap", label = T)
res0.9_metadata <- boule_subsetted@meta.data
#trim the metadata containing the boule cluster cells
res0.9_metadata <- res0.9_metadata[,-4]

#writing the metadata file into the directory of choice.
write.table(res0.9_metadata, "/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/res0.9_metadata.txt", sep="\t")

#read in custom metadata after merging the metadatas together using a python script
custom_metadata <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Integrated_Regeneration/JK_embed_integrated_regeneration/custom_metadata.txt", sep = "\t")

#take the column with the custom clustering information and store it in an object called "res"
res <- custom_metadata[,4]

#take the object "res" and make it a new column in the integrated dataset's metadata
int_reg_neoblast@meta.data <- cbind(int_reg_neoblast@meta.data, res)

#set the active identity for the int dataset to this "res" column so it shows up in the umap and other commands
int_reg_neoblast <- SetIdent(int_reg_neoblast, value = int_reg_neoblast@meta.data$res)
int_reg_neoblast@active.ident

#plot umap with the "res" clustering column by specifying the group.by parameter.
DimPlot(object = int_reg_neoblast, reduction = "umap",group.by = "res", label = T)




