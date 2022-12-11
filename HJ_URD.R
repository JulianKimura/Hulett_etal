###########PREP WORK FOR RUNNING URD ON A SINGLE DATASET

#####First, try looking in to pseudotime analysis using scanpy. 
#need to obtain the embedding data through here:

embeddings <- juv.data@reductions$umap@cell.embeddings
write.table(embeddings, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/Trajectory/HJ_cell_embedings.txt", sep = '\t')

HJ_matrix <- juv.data@assays$RNA@counts
write.table(HJ_matrix, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/Trajectory/HJ_matrix.txt", sep = '\t')

HJ_metadata <- juv.data@meta.data
write.table(HJ_metadata, "~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/Trajectory/HJ_metadata.txt", sep = '\t')




########now that we have an idea on where the clusters are in pseudotime we can subset the original seurat object into different parts to be processed as different stages in URD.
library(Seurat)

#first load the R data for the juvenile seurat clustering. 
#below are the parameters that were used to create the UMAP

#run the dimplot to show the UMAP clustering
DimPlot(object = juv.data, reduction = "umap", label = T, group.by = 'res')
juv.data <- SetIdent(juv.data, value = juv.data@meta.data$res)
juv.data@active.ident

#subset the dataset into different cell populations (terminal cell types, neoblasts, and intermediate)
term_cell<- subset(juv.data, idents =c("0","19", "5", "14", "10", "1", "12", "2", "6",  "8", "16", "15", "20", "17", "18", "21"))

DimPlot(object = term_cell, reduction = "umap", label = T)


neoblast<- subset(juv.data, idents =c("3"))

DimPlot(object = neoblast, reduction = "umap", label = T)


int_cell<- subset(juv.data, idents =c("4","9","13" ,"11","7"))

DimPlot(object = int_cell, reduction = "umap", label = T)

#extract cell population specific matrices
term_cell.matrix <- as.matrix(term_cell@assays$RNA@data)
neoblast.matrix <- as.matrix(neoblast@assays$RNA@data)
int_cell.matrix <- as.matrix(int_cell@assays$RNA@data)

#extract cell population specific metadata
term_cell.metadata <- term_cell@meta.data
term_cell.metadata <- term_cell.metadata[,-5]
term_cell.metadata$orig.ident <- "terminal"
#names(term_cell.metadata)[names(term_cell.metadata) == "RNA_snn_res.0.5"] <- "res"

neoblast.metadata <- neoblast@meta.data
neoblast.metadata <- neoblast.metadata[,-5]
neoblast.metadata$orig.ident <- "neoblast"
#names(neoblast.metadata)[names(neoblast.metadata) == "RNA_snn_res.0.5"] <- "res"


int_cell.metadata <- int_cell@meta.data
int_cell.metadata <- int_cell.metadata[,-5]
int_cell.metadata$orig.ident <- "intermediate"
#names(int_cell.metadata)[names(int_cell.metadata) == "RNA_snn_res.0.5"] <- "res"


######combining matrices
list_cells <- list(neoblast.matrix,term_cell.matrix,int_cell.matrix)
combined_matrix <- dsCombineDGE(list_cells)

combined_metadata <- do.call("rbind", list(neoblast.metadata,int_cell.metadata,term_cell.metadata))


###########RUNNING URD
library(Seurat)
library(magrittr)
library(dplyr)
require(scales)
library(URD)
library(cowplot)


# Create an URD object, which will filter the data, then normalize and log-transform it.
embryo <- createURD(count.data = combined_matrix, meta = combined_metadata, min.cells=5, min.counts=200)

# Copy stage from @meta to @group.ids 
#this takes in the stage.nice columns in the metadata, and placing it into a new slot in the axial data object called stage within group.ids. 
embryo@group.ids$stage <- as.character(embryo@meta[rownames(embryo@group.ids),"orig.ident"])
embryo@group.ids
# Get variable genes for each group of 3 stages
# (Normally would do this for each stage, but there are not very many cells in this subset of the data)
# diffCV.cutoff can be varied to include more or fewer genes.

#now you are taking the stage names that you had in your group.ids slot in the URD object, and creating a new object called stages. This takes all of the unique stages and then creates
stages <- sort(unique(embryo@group.ids$stage))
stages

#you are taking the variable genes per stage here. 
var.by.stage <- lapply(seq(1,3,1), function(n) {findVariableGenes(embryo, cells.fit=cellsInCluster(embryo, "stage", stages), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, main.use=paste0("Stages ", stages[n-2], " to ", stages[n]), do.plot=T)})


save.image("~/Documents/Lab/Single_Cell_Seq/juv_brov2_pca.RData")


# Combine the results from each group of stages into a single list of variable genes and load into the URD object
var.genes <- sort(unique(unlist(var.by.stage)))
embryo@var.genes <- var.genes

# Calculate PCA and consider those PCs that with standard deviation 2x expected by noise as significant
embryo <- calcPCA(embryo, mp.factor = 2)

pcSDPlot(embryo)


# Calculate tSNE
set.seed(19)
embryo <- calcTsne(object = embryo)
plotDim(embryo, "orig.ident", plot.title = "tSNE: Stage")

# In this case, knn=100 (larger than sqrt(n.cells)) works well because there are not many cell types.
# Sigma 16 is slightly smaller than the sigma auto-determined by using NULL parameter.
embryo <- calcDM(embryo, knn = 91, sigma=NULL)


plotDimArray(embryo, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map(Sigma NULL 91NNs): Stage", label="stage", plot.title="", legend=F)


plotDim(embryo, "orig.ident", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")


#calculate pseudotime
# Here we use all cells from the first stage as the root
root.cells <- cellsInCluster(embryo, "stage", "neoblast")


########THIS PART OF ANALYSIS NEEDS TO BE DONE ON CLUSTER
# Then we run 'flood' simulations
embryo.floods <- floodPseudotime(embryo, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

# The we process the simulations into a pseudotime
embryo <- floodPseudotimeProcess(embryo, embryo.floods, floods.name="pseudotime")


pseudotimePlotStabilityOverall(embryo)
plotDim(embryo, "pseudotime")

pdf(file="/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/Embryo_juv_v2.0/pseudotime.pdf",height=4.1,width=9.8)
p5<-plotDim(embryo, "pseudotime")
plotgrid(p5)
dev.off()

pdf(file="/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/Embryo_juv_v2.0/pseudotime_by_stage.pdf",height=4.1,width=9.8)
p6<-plotDists(embryo, "pseudotime", "stage", plot.title="Pseudotime by stage")
plotgrid(p6)
dev.off()

plotDists(embryo, "pseudotime", "res", plot.title="Pseudotime by cluster", label = 'res')

plotDim(embryo, "pseudotime")

# Create a subsetted object of just those cells from the final stage
embryo.juv <- urdSubset(embryo, cells.keep=cellsInCluster(embryo, "stage", "terminal"))

# Use the variable genes that were calculated only on the final group of stages (which
# contain the last stage).
embryo.juv@var.genes <- var.by.stage[[3]]

# Calculate PCA and tSNE
embryo.juv <- calcPCA(embryo.juv, mp.factor = 1.5)
pdf(file="/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/Embryo_juv_v2.0/pcSDPlot.pdf",height=4.1,width=9.8)
p7<-pcSDPlot(embryo.juv)
plotgrid(p7)
dev.off()

set.seed(20)
embryo.juv <- calcTsne(embryo.juv, which.dims = 1:10, perplexity = 50)

pdf(file="/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/Embryo_juv_v2.0/seurat_int_clustering.pdf",height=4.1,width=9.8)
p8<-plotDim(embryo.juv, "res", plot.title = "seurat_J_clustering", point.size=3)
plotgrid(p8)
dev.off()

#biased random walks
# Copy cluster identities from axial.6somite object to a new clustering ("tip.clusters") in the full axial object.
embryo@group.ids[rownames(embryo.juv@group.ids), "tip.clusters"] <- embryo.juv@meta$`res`

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
embryo.ptlogistic <- pseudotimeDetermineLogistic(embryo, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)

# Bias the transition matrix acording to pseudotime
embryo.biased.tm <-as.matrix(pseudotimeWeightTransitionMatrix(embryo, "pseudotime", logistic.params=embryo.ptlogistic))

# Simulate the biased random walks from each tip
embryo.walks <- simulateRandomWalksFromTips(embryo, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = embryo.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)

# Process the biased random walks into visitation frequencies
embryo <- processRandomWalksFromTips(embryo, embryo.walks, verbose = F)

# Load the cells used for each tip into the URD object
embryo.tree <- loadTipCells(embryo, "tip.clusters")

# Build the tree
embryo.tree <- buildTree(embryo.tree, pseudotime = "pseudotime", tips.use= c("0","19", "5", "14", "10", "1", "12", "2", "6",  "8", "16", "15", "20", "17", "18", "21"), divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

plotTree(embryo.tree, "segment", label.segments = T)
plotTree(embryo.tree, "stage", label.segments = T)

pdf(file="/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/Embryo_juv_v2.0/embryo_tree.pdf",height=9.8,width=9.8)
p9<-plotTree(embryo.tree, "stage", title="Developmental Stages")
plot_grid(p9)
dev.off()


save.image("/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/URD_v2.0_Tree.RData")


embryo.tree <- nameSegments(embryo.tree, segments=c("0","1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"), segment.names = c("0","1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

###############

plotTree(embryo.tree, "segment", title="URD Tree Segment", label.segments = TRUE)
plotTree(embryo.tree, "stage", title="URD Tree Stages")



embryo.muscle <- urdSubset(embryo, cells.keep=cellsInCluster(embryo, "segment", "25"))

embryo.juv <- urdSubset(embryo, cells.keep=cellsInCluster(embryo, "stage", "terminal"))



#Warning message:
#In pseudotimeBreakpointByStretch(div.pseudotime, segment.1, segment.2,  :
#No obvious breakpoint between 15 and 27 . Longest stretch of difference is #upstream of longest stretch of non-different.


#########force directed layout
embryo.tree <- treeForceDirectedLayout(embryo.tree, num.nn=NULL, cut.unconnected.segments=1, verbose=T)
plotTreeForce(embryo.tree, "segment")

embryo.tree <- treeForceRotateCoords(embryo.tree, seg="1", angle = -3.3, axis="z", around.cell = 10, throw.out.cells=1000, pseudotime="pseudotime")



plotTreeForceStore3DView(embryo.tree, view1)
plotTreeForce(embryo.tree, "segment")
plotTreeForce(embryo.tree, "stage")

plotTreeForce(embryo.tree, "98046269|ryr1")

plotTreeForce2D(embryo.tree, label = "segment", label.type = "search", show.points = T, point.alpha = 1, point.size = 1,show.neighbors = T, neighbors.max = 10000, colors = NULL)

plotTreeForce(embryo.tree, "98103427|abcb6-6")


plotTree(embryo.tree, "98017127|pax5", title="98017127|pax5")

plotTree(embryo.tree, "98008245-myof", title="98008245-myof")
plotTree(embryo.tree, "98008634-foxf1", title="98008634-foxf1")
plotTree(embryo.tree, "98005071-vax1", title="98005071-vax1")
plotTree(embryo.tree, "98003853-dlx1", title="98003853-dlx1")
plotTree(embryo.tree, "98007047-ikzf2", title="98007047-ikzf2")
plotTree(embryo.tree, "98033897-six1", title="98033897-six1")
plotTree(embryo.tree, "98054882-foxj1", title="98054882-foxj1")
plotTree(embryo.tree, "98045757-nkx24-2", title="98045757-nkx24-2")
plotTree(embryo.tree, "98037905-sox4", title="98037905-sox4")


# Determine tips to run DE for


###markersaucpr is a way to take two cell pops and ask diff genes. 
###aucprtestalong tree calls this multiple times, and finds markers.

##this labels the segments!!!
plotTree(embryo.tree, "segment", label.segments = T)

#decide what tips to run DE for
tips.to.run <- c("8")
#genes.use <- NULL # Calculate for all genes
genes.use <- embryo.tree@var.genes
# Calculate the markers of each other population.
gene.markers <- list()
for (tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[tipn]
  print(paste0(Sys.time(), ": ", tip))
  markers <- aucprTestAlongTree(embryo.tree, pseudotime="pseudotime", tips=tip, log.effect.size=0.4, auc.factor = 1.25, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use=genes.use, root="35", only.return.global=F, must.beat.sibs=0.1, report.debug=T, segs.to.skip = NULL)
  saveRDS(markers, file=paste0("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/trajectory/Neural/tip_8_root_35", tip, "tip25_root35.rds"))
  gene.markers[[tip]] <- markers
}
write.table(gene.markers[[1]]$diff.exp, file="/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/trajectory/Neural/tip_8_root_35.txt", sep="\t")



#this gives you top genes that have higher expression in relation to the other branches. 
head(gene.markers[[1]]$marker.chain)

#this gives you the differentially expressed genes
head(gene.markers[[1]]$diff.exp)
write.table(gene.markers[[1]]$diff.exp, file="/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/INTEGRATED_FILTERED/tip11_root36.txt", sep="\t")

#mustbeatsibs: if looking for marker of 0, it compares to all of its siblings, what percentage of its signlings it must beat. great for polytomys. 

#calculate markers, loop through, then save them sort through many different parameters.


#gene.markers[[1]]$marker.chain
#gene.markers[[1]]$stats


#The primordial germ cells (PGCs) and enveloping layer cells (EVL) are distinct from the beginning of our tree. Thus, they are not divided up into small segments, and our differential expression test along the tree performed poorly for them. So, instead, we divided cells into five groups, based on their developmental stage, and performed pairwise comparisons at each stage using the precision-recall approach, and kept those genes that were markers of at least 3 groups.

####For pulling out DE based on stage in segments
# Define cell populations
neoblast.cells <- cellsInCluster(embryo.tree, "segment", "2")
aneoblast2.cells <- cellsInCluster(embryo.tree, "segment", "12")
root.cells <- cellsInCluster(object, "segment", segChildrenAll(embryo.tree, "27", include.self=T))
# Copy STAGE to group.ids for the TestByFactor function
embryo.tree@group.ids$STAGE <- embryo.tree@meta[rownames(embryo.tree@group.ids), "STAGE"]
# Define stage groups
groups <- list(
  c("E80", "ZFOBLONG"),
  c("ZFDOME", "ZF30"),
  c("ZF50", "ZFS", "ZF60"),
  c("ZF75", "ZF90"),
  c("ZFB", "ZF3S", "ZF6S")
)
# Calculate markers
evl.markers.bystage <- aucprTestByFactor(object, cells.1=evl.cells, cells.2=list(pgc.cells, blastoderm.cells),
                                         label="STAGE", groups=groups, 
                                         log.effect.size=0.5, auc.factor=1, min.auc.thresh=0.1, max.auc.thresh=Inf,
                                         frac.must.express=0.1, frac.min.diff=0, genes.use=genes.use, min.groups.to.mark=3, report.debug=T)
pgc.markers.bystage <- aucprTestByFactor(object, cells.1=pgc.cells, cells.2=list(evl.cells, blastoderm.cells),
                                         label="STAGE", groups=groups, 
                                         log.effect.size=0.5, auc.factor=1, min.auc.thresh=0.1, max.auc.thresh=Inf,
                                         frac.must.express=0.1, frac.min.diff=0, genes.use=genes.use, min.groups.to.mark=3, report.debug=T)
# Save them
saveRDS(evl.markers.bystage, "cascades/aucpr/EVL/Periderm.rds")
saveRDS(pgc.markers.bystage, "cascades/aucpr/Primordial Germ Cells.rds")



#force directed layout
# Generate the force-directed layout
early.tree <- treeForceDirectedLayout(early.tree, num.nn=100, cut.unconnected.segments=2, verbose=T)

plotTreeForce(axial.tree, "GSC", title = "GSC", title.cex = 2, title.line=2.5)
plotTreeForce(axial.tree, "HE1A", title = "HE1A", title.cex=2, title.line=2.5)
plotTreeForce(axial.tree, "COL8A1A", title="COL8A1A", title.cex=2, title.line=2.5)

pdf(file="/n/srivastava_lab/Julian_Kimura/URD_EMBRYO/INTEGRATED_FILTERED/logistic.pdf",height=9.8,width=9.8)



#####genes that are good markers for the neoblast specification. 
plotTree(embryo.tree, "98007047-ikzf2", title="ikzf2")
plotTree(embryo.tree, "98027755-foxa1", title="foxa1")
plotTree(embryo.tree, "98003853-dlx1", title="dlx1")
plotTree(embryo.tree, "98008634-foxf1", title="foxf1")
plotTree(embryo.tree, "98054882-foxj1", title="foxj1")
plotTree(embryo.tree, "98005071-vax1", title="vax1")
plotTree(embryo.tree, "98033897-six1", title="six1")


plotTree(embryo.tree, "98027755-foxa1", title="98027755-foxa1")







#######additional work to make the HJ data UMAP cleaner

HJ <- JackStraw(HJ, num.replicate = 100)
HJ <- ScoreJackStraw(HJ, dims = 1:20)
JackStrawPlot(HJ, dims = 1:20)


HJ <- FindNeighbors(object = HJ, dims = 1:15)
HJ <- FindClusters(HJ, resolution = 0.5, print.output = 0, save.SNN = T)
HJ <-  RunUMAP(HJ, reduction.use = "pca", dims= 1:20, min.dist = 0.6, n.neighbors = 50)

DimPlot(object = HJ, reduction = "umap",label = T)




########generating HJ pseudobulk matrices.


#########subsetting the tree based off of particular lineages
plotTree(embryo.tree, "segment", title="URD Tree Segment", label.segments = TRUE)
plotTree(embryo.tree, "stage", title="URD Tree Segment", label.segments = F)


hj_int_subset <- urdSubset(embryo.tree, cells.keep=cellsInCluster(embryo.tree, clustering = "segment", cluster = c("35")))

#intermediate_subset <- urdSubset(intermediate_subset, cells.keep=cellsInCluster(embryo.tree, clustering = "stage", cluster = c("E35", "E50", "E65", "E80", "E95", "E110", "E125", "E145")))

plotTree(hj_int_subset, "stage", title="URD Tree Stages", plot.tree =T)


hj_int_clade_matrix <- as.matrix(hj_int_subset@count.data)

write.table(hj_int_clade_matrix, "~/Documents/Lab/embryo_singlecell_paper/hj_int_clade_matrix.txt", sep = "\t")



#########create a new subsetted tree that is missing the juvenile dataset, that way I can then perfom differential expression analysis for specific stages and clusters. 


embryo.tree.only <- urdSubset(embryo.tree, cells.keep=cellsInCluster(embryo.tree, clustering = "stage", cluster = c("E35", "E50", "E65", "E80", "E95", "E110", "E125", "E145")))

plotTree(embryo.tree.only, "stage", title="URD Tree Stages", plot.tree =T, label.segments = T)

#now do differential expression for each of the clades that I created pseudobulk matrices for. 


plotTree(embryo.tree, "stage", title="URD Tree Segment", label.segments = TRUE)

#decide what tips to run DE for
tips.to.run <- c("0")
#genes.use <- NULL # Calculate for all genes
genes.use <- embryo.tree@var.genes
# Calculate the markers of each other population.
gene.markers <- list()
for (tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[tipn]
  print(paste0(Sys.time(), ": ", tip))
  markers <- aucprTestAlongTree(embryo.tree, pseudotime="pseudotime", tips=tip, log.effect.size=0.4, auc.factor = 1.25, max.auc.threshold = 0.1, frac.must.express = 0.1, frac.min.diff = 0, genes.use=genes.use, root="35", only.return.global=F, must.beat.sibs=0.1, report.debug=T, segs.to.skip = NULL)
  saveRDS(markers, file=paste0("/Users/JulianKimura/Documents/Lab/embryo_singlecell_paper/hj_clade_matrices/", tip, "neural_root35.rds"))
  gene.markers[[tip]] <- markers
}
write.table(gene.markers[[1]]$diff.exp, file="/Users/JulianKimura/Documents/Lab/embryo_singlecell_paper/hj_gut_root35.txt", sep="\t")



plotTree(embryo.tree, "98030438-hxb1", title = "")



