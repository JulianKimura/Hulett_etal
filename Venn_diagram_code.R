
library(VennDiagram)

setwd("~/Documents/Lab/Single_Cell_Seq/Post_Embryonic/RH_For_JK_URD_Data2/JK_embedding/Figures/")

#big question for the venn diagram: do we count the repeated cell clusters as independent ones? for example, we have two clusters that corresponds to mesenchymal based on HJ clustering. do we keep this terminology? 

HJ <- paste(1:10)
HJ2 <- paste(c("unknown v", "unknown vii", "unknown vi", "unknown ii", "unknown i", "unknown viii", "unknown ix"))
#HJ2 <- paste(rep("HJ", 7))
HJ <- c(HJ, HJ2)
rm(HJ2)

LJ <- paste(1:10)
LJ2 <- paste(c("unknown i", "unknown viii", "unknown ii", "unknown iv", "unknown ix"))
#LJ2 <- paste(rep("LJ", 5))
LJ <- c(LJ, LJ2)
rm(LJ2)


EA <- paste(1:10)
EA2 <- paste(c("unknown v", "unknown ii", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"))
#EA2 <- paste(rep("EA", 13))
EA <- c(EA, EA2)
rm(EA2)

A <- paste(1:10)
A2 <- paste(c("unknown ix", "11", "12", "13", "14", "15", "16", "22", "23", "24", "25", "26", "27"))
#A2 <- paste(rep("A", 13))
A <- c(A, A2)
rm(A2)

myCol <- brewer.pal(4, "Pastel2")

venn.diagram(x = list(HJ, LJ, EA, A), category.names = c("Hatchling Juvenile", "Late Juvenile", "Early Adult", "Adult"), filename = "test_cell_type_venn_diagram.png", output = T, imagetype="png" ,height = 2000 , width = 3000 , resolution = 300,compression = "lzw",fill = myCol)



test <- list(A,EA,LJ,HJ)

ggVennDiagram(test, label_alpha = 0)


ggVennDiagram(
  test, label_alpha = 0,
  category.names = c("Adult","Early Adult","Late Juvenile", "Hatchling Juvenile")
) +
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")

