####This Rscript should only be run after you have generated the file titled hmi_cell_cluster_embryonic.txt


# Load library
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
# Load the matrix you  just generated using the bash commands.
matrix_data <- read.table("/Users/JulianKimura/Documents/Lab/embryo_singlecell_paper/clade_pearsons/adding_pipes/hmi_cell_cluster_embryonic.txt", header=T)
# Calculate Pearson's
matrix_cor <- cor(matrix_data, use ="all.obs" ,method ="pearson")
# Setup color palette
my_palette <- colorRampPalette(c("royal blue", "white", "red"))(n = 128)
# Draw heatmap
pheatmap(matrix_cor, fontsize_row=5, fontsize_col=5, cluster_cols=T, cluster_rows=T,
         col=my_palette, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="non
 e", cexCol=1, cexRow=1, cellwidth=5, cellheight=5)

