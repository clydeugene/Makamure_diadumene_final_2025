# .libPaths(c("/home/samir/R/x86_64-pc-linux-gnu-library/4.0/", .libPaths()))
# Load required packages
library("ComplexHeatmap")
library("circlize")

# Load data
args <- as.character(commandArgs(trailingOnly=T))
NRC_mat <- read.delim(args[1],header=T,row.names=1) # Header may not be present
Log2_comp <- NRC_mat

head(Log2_comp[,1])

Log2_mat <- as.matrix(read.table(args[2],header=T,row.names=1))
Hclust_reps <- Log2_mat

# Read heatmap name from arguments
heatmap_name <- args[3]

# Match row order
Hclust_reps <- Hclust_reps[rownames(Log2_comp),]

# Select the appropriate column input for the Log2 heatmap
Hclust_rep_RT_1 <- Hclust_reps[,grep("Room_1", colnames(Hclust_reps))]
Hclust_rep_RT_2 <- Hclust_reps[,grep("Room_2", colnames(Hclust_reps))]
Hclust_rep_RT_3 <- Hclust_reps[,grep("Room_3", colnames(Hclust_reps))]
Hclust_rep_RT_4 <- Hclust_reps[,grep("Room_4", colnames(Hclust_reps))]
Hclust_rep_RT_5 <- Hclust_reps[,grep("Room_5", colnames(Hclust_reps))]

Hclust_rep_16_1 <- Hclust_reps[,grep("Cold_1", colnames(Hclust_reps))]
Hclust_rep_16_2 <- Hclust_reps[,grep("Cold_2", colnames(Hclust_reps))]
Hclust_rep_16_3 <- Hclust_reps[,grep("Cold_3", colnames(Hclust_reps))]
Hclust_rep_16_4 <- Hclust_reps[,grep("Cold_4", colnames(Hclust_reps))]
Hclust_rep_16_5 <- Hclust_reps[,grep("Cold_5", colnames(Hclust_reps))]

head(Log2_comp[,1])

# Create heatmaps for Hclust TPM values
heatmap_Hclust_rep_RT_1 <- Heatmap(Hclust_rep_RT_1, name = "RT_1", col = colorRamp2(c(1, 12), c( "white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_2 <- Heatmap(Hclust_rep_RT_2, name = "RT_2", col = colorRamp2(c(1, 12), c( "white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_3 <- Heatmap(Hclust_rep_RT_3, name = "RT_3", col = colorRamp2(c(1, 12), c( "white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_4 <- Heatmap(Hclust_rep_RT_4, name = "RT_4", col = colorRamp2(c(1, 12), c( "white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_5 <- Heatmap(Hclust_rep_RT_5, name = "RT_5", col = colorRamp2(c(1, 12), c( "white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)

heatmap_Hclust_rep_16_1 <- Heatmap(Hclust_rep_16_1, name = "16_1", col = colorRamp2(c(1, 12), c("white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_2 <- Heatmap(Hclust_rep_16_2, name = "16_2", col = colorRamp2(c(1, 12), c("white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_3 <- Heatmap(Hclust_rep_16_3, name = "16_3", col = colorRamp2(c(1, 12), c("white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_4 <- Heatmap(Hclust_rep_16_4, name = "16_4", col = colorRamp2(c(1, 12), c("white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_5 <- Heatmap(Hclust_rep_16_5, name = "16_5", col = colorRamp2(c(1, 12), c("white", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)

head(Hclust_rep_RT_1)
head(Log2_comp[,1])

col_fun <- colorRamp2(c(1, 12), c("white", "black"))

# Create heatmaps for Log2 values
heatmap_Log2_comp_16_vs_RT_2 <- Heatmap(Log2_comp[,1], name = "Log2Fold", col = colorRamp2(c(-3, 0, 3), c("blue","white", "red")), show_row_names = FALSE, cluster_rows = FALSE)

# Draw the combined heatmap
png(paste(heatmap_name, ".png", sep = ""), height = 2000, width = 2000, res = 300)
draw(heatmap_Log2_comp_16_vs_RT_2 +
    heatmap_Hclust_rep_RT_1 +
    heatmap_Hclust_rep_RT_2 +
    heatmap_Hclust_rep_RT_3 +
    heatmap_Hclust_rep_RT_4 +
    heatmap_Hclust_rep_RT_5 +
    heatmap_Hclust_rep_16_1 +
    heatmap_Hclust_rep_16_2 +
    heatmap_Hclust_rep_16_3 +
    heatmap_Hclust_rep_16_4 +
    heatmap_Hclust_rep_16_5
)
dev.off()


# Make colorful version of the heatmap
heatmap_Hclust_rep_RT_1 <- Heatmap(Hclust_rep_RT_1, name = "RT_1", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c( "white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_2 <- Heatmap(Hclust_rep_RT_2, name = "RT_2", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c( "white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_3 <- Heatmap(Hclust_rep_RT_3, name = "RT_3", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c( "white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_4 <- Heatmap(Hclust_rep_RT_4, name = "RT_4", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c( "white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_RT_5 <- Heatmap(Hclust_rep_RT_5, name = "RT_5", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c( "white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)

heatmap_Hclust_rep_16_1 <- Heatmap(Hclust_rep_16_1, name = "16_1", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c("white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_2 <- Heatmap(Hclust_rep_16_2, name = "16_2", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c("white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_3 <- Heatmap(Hclust_rep_16_3, name = "16_3", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c("white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_4 <- Heatmap(Hclust_rep_16_4, name = "16_4", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c("white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Hclust_rep_16_5 <- Heatmap(Hclust_rep_16_5, name = "16_5", col = colorRamp2(c(1, 3.75, 6.5, 9.25, 12), c("white", "green", "purple", "orange", "black")), heatmap_legend_param = list(at = c(1, 3.75, 6.5, 9.25, 12)), show_row_names = FALSE, cluster_rows = FALSE)

head(Hclust_rep_RT_1)
head(Log2_comp[,1])

col_fun <- colorRamp2(c(1, 12), c("white", "black"))

# Create heatmaps for Log2 values
heatmap_Log2_comp_16_vs_RT_2 <- Heatmap(Log2_comp[,1], name = "Log2Fold", col = colorRamp2(c(-3, 0, 3), c("blue","white", "red")), show_row_names = FALSE, cluster_rows = FALSE)

# Draw the combined heatmap
png(paste(heatmap_name, "colorful", ".png", sep = ""), height = 2000, width = 2000, res = 300)
draw(heatmap_Log2_comp_16_vs_RT_2 +
    heatmap_Hclust_rep_RT_1 +
    heatmap_Hclust_rep_RT_2 +
    heatmap_Hclust_rep_RT_3 +
    heatmap_Hclust_rep_RT_4 +
    heatmap_Hclust_rep_RT_5 +
    heatmap_Hclust_rep_16_1 +
    heatmap_Hclust_rep_16_2 +
    heatmap_Hclust_rep_16_3 +
    heatmap_Hclust_rep_16_4 +
    heatmap_Hclust_rep_16_5
)
dev.off()
