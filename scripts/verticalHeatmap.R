library("ComplexHeatmap")
library("circlize")
library("pheatmap")
library("fastcluster")

args=as.character(commandArgs(trailingOnly=T))
M=read.delim(args[1],header=T,row.names=1)
Compare=as.character(args[2])

head(M[,1])

# Count number of DE genes/isoforms
length(M[,1])

N=M
head(N)
ordr=hclust(dist(N))$order

H=Heatmap(M[ordr,], name = paste("TPM", sep=""), col = colorRamp2(c(-1.75, -.2, .2, 1.75), c("blue", "white", "white", "red")), show_row_names = FALSE, cluster_rows = FALSE, show_heatmap_legend = TRUE, cluster_columns = FALSE)

png(paste(Compare, "_heatmap_vertical.png", sep=""), height = 2000, width = 2000, res = 300, type = "cairo")
draw(H)
dev.off()


write.table(M[ordr,], file="Hclust_TPMTable_2.txt", sep="\t", quote=FALSE)
