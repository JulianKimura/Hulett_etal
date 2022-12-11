#YJLs code for making dotplots
# ggplot2 point plot (pp) with specified circle size
library(ggplot2)

# GO_BP
# ggplot2 point plot (pp) with specified circle size (5x5.5 inches)
data <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/GO_dotplots/dotplot_matrix.txt", sep = '\t', header = T)
pp <- ggplot(data, aes(Cluster,GO.Terms))
pp <- pp + geom_point(aes(size = Counts, colour= -log(FDR)), alpha=1.0) + scale_size(range = c(6, 12))
pp <- pp + scale_colour_gradientn(colours = c("royal blue","light grey","red"))
pp <- pp + scale_x_discrete(limits=c("C0", "C1", "C2", "C3", "C5"), breaks = c("C0", "C1", "C2", "C3", "C5"))
pp <- pp + scale_y_discrete(limits=c("lipid metabolic process", "muscle system process", "synaptic vesicle transport", "mRNA metabolic process", "cilium organization"))
pp <- pp + theme(panel.grid.major = element_line(colour = "light grey",size=0.2), panel.grid.minor = element_blank(), axis.ticks=element_line(colour = "black"),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 axis.text=element_text(colour = "black", size=6), axis.title=element_text(size=8))
pp <- pp + theme(panel.border=element_blank(), axis.line=element_line())

pp




data <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/HJ/filtered/GO_dotplots/dotplot_matrix.txt", sep = '\t', header = T)
pp <- ggplot(data, aes(Cluster,GO.Terms))
pp <- pp + geom_point(aes(size = Counts, colour= -log(FDR)), alpha=1.0) + scale_size(range = c(6, 12))
pp <- pp + scale_colour_gradientn(colours = c("royal blue","light grey","red"))
pp <- pp + scale_x_discrete(limits=data$Cluster, breaks = data$Cluster)
pp <- pp + scale_y_discrete(limits=data$GO.Terms)
pp <- pp + theme(panel.grid.major = element_line(colour = "light grey",size=0.2), panel.grid.minor = element_blank(), axis.ticks=element_line(colour = "black"),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 axis.text=element_text(colour = "black", size=6), axis.title=element_text(size=8))
pp <- pp + theme(panel.border=element_blank(), axis.line=element_line())

pp
