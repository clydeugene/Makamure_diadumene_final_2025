# Load required packages
library("ComplexHeatmap")
library("circlize")


# Load data
args <- as.character(commandArgs(trailingOnly = TRUE))

log2_comp <- read.delim(args[1], header = FALSE, row.names = 1)

log2_mat <- as.matrix(
  read.table(args[2], header = TRUE, row.names = 1, check.names = FALSE)
)

heatmap_name <- args[3]


hclust_reps <- log2_mat

# Match row order
hclust_reps <- hclust_reps[rownames(log2_comp), ]

# Select the appropriate column input for the log2 heatmap
hclust_rep_rt_1 <- hclust_reps[,grep("rt_rep_1", colnames(hclust_reps))]
hclust_rep_rt_2 <- hclust_reps[,grep("rt_rep_2", colnames(hclust_reps))]
hclust_rep_rt_3 <- hclust_reps[,grep("rt_rep_3", colnames(hclust_reps))]
hclust_rep_rt_4 <- hclust_reps[,grep("rt_rep_4", colnames(hclust_reps))]
hclust_rep_rt_5 <- hclust_reps[,grep("rt_rep_5", colnames(hclust_reps))]

hclust_rep_16_1 <- hclust_reps[,grep("16_rep_1", colnames(hclust_reps))]
hclust_rep_16_2 <- hclust_reps[,grep("16_rep_2", colnames(hclust_reps))]
hclust_rep_16_3 <- hclust_reps[,grep("16_rep_3", colnames(hclust_reps))]
hclust_rep_16_4 <- hclust_reps[,grep("16_rep_4", colnames(hclust_reps))]
hclust_rep_16_5 <- hclust_reps[,grep("16_rep_5", colnames(hclust_reps))]

head(log2_comp[, 1])

col_fun <- colorRamp2(c(1, 12), c("white", "black"))
heat_leg <- list(at = c(1, 3.75, 6.5, 9.25,  12))

# Create heatmaps for hclust TPM values
heatmap_hclust_rep_rt_1 <- Heatmap(hclust_rep_rt_1, name = "rt_rep_1", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_2 <- Heatmap(hclust_rep_rt_2, name = "rt_rep_2", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_3 <- Heatmap(hclust_rep_rt_3, name = "rt_rep_3", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_4 <- Heatmap(hclust_rep_rt_4, name = "rt_rep_4", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_5 <- Heatmap(hclust_rep_rt_5, name = "rt_rep_5", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)

heatmap_hclust_rep_16_1 <- Heatmap(hclust_rep_16_1, name = "16_rep_1", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_2 <- Heatmap(hclust_rep_16_2, name = "16_rep_2", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_3 <- Heatmap(hclust_rep_16_3, name = "16_rep_3", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_4 <- Heatmap(hclust_rep_16_4, name = "16_rep_4", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_5 <- Heatmap(hclust_rep_16_5, name = "16_rep_5", col = col_fun, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)

head(hclust_rep_rt_1)
head(log2_comp[, 1])

# Create heatmaps for log2 values
heatmap_log2_comp_16_vs_rt_2 <- Heatmap(log2_comp[,1], name = "log2Fold", col = colorRamp2(c(-3, 0, 3), c("blue","white", "red")), show_row_names = FALSE, cluster_rows = FALSE)

# Draw the combined heatmap
png(paste(heatmap_name, ".png", sep = ""), height = 2000, width = 2000, res = 300)
draw(heatmap_log2_comp_16_vs_rt_2 +
    heatmap_hclust_rep_rt_1 +
    heatmap_hclust_rep_rt_2 +
    heatmap_hclust_rep_rt_3 +
    heatmap_hclust_rep_rt_4 +
    heatmap_hclust_rep_rt_5 +
    heatmap_hclust_rep_16_1 +
    heatmap_hclust_rep_16_2 +
    heatmap_hclust_rep_16_3 +
    heatmap_hclust_rep_16_4 +
    heatmap_hclust_rep_16_5
)
dev.off()

# Save the heatmaps as a pdf
pdf(paste(heatmap_name, ".pdf", sep = ""), height = 8, width = 8)
draw(heatmap_log2_comp_16_vs_rt_2 +
    heatmap_hclust_rep_rt_1 +
    heatmap_hclust_rep_rt_2 +
    heatmap_hclust_rep_rt_3 +
    heatmap_hclust_rep_rt_4 +
    heatmap_hclust_rep_rt_5 +
    heatmap_hclust_rep_16_1 +
    heatmap_hclust_rep_16_2 +
    heatmap_hclust_rep_16_3 +
    heatmap_hclust_rep_16_4 +
    heatmap_hclust_rep_16_5
)
dev.off()

################################################
############ Make colorful heatmap #############
################################################


# Define a new color scheme for the heatmap
col_fun_2 <- colorRamp2(
  c(1, 3.75, 6.5, 9.25, 12),
  c("white", "green", "purple", "orange", "black")
)

# Create individual heatmaps for each replicate
heatmap_hclust_rep_rt_1 <- Heatmap(hclust_rep_rt_1, name = "rt_rep_1", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_2 <- Heatmap(hclust_rep_rt_2, name = "rt_rep_2", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_3 <- Heatmap(hclust_rep_rt_3, name = "rt_rep_3", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_4 <- Heatmap(hclust_rep_rt_4, name = "rt_rep_4", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_rt_5 <- Heatmap(hclust_rep_rt_5, name = "rt_rep_5", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)

heatmap_hclust_rep_16_1 <- Heatmap(hclust_rep_16_1, name = "16_rep_1", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_2 <- Heatmap(hclust_rep_16_2, name = "16_rep_2", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_3 <- Heatmap(hclust_rep_16_3, name = "16_rep_3", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_4 <- Heatmap(hclust_rep_16_4, name = "16_rep_4", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)
heatmap_hclust_rep_16_5 <- Heatmap(hclust_rep_16_5, name = "16_rep_5", col = col_fun_2, heatmap_legend_param = heat_leg, show_row_names = FALSE, cluster_rows = FALSE)


# Create heatmap for log2 values
heatmap_log2_comp_16_vs_rt_2 <- Heatmap(log2_comp[,1], name = "log2Fold", col = colorRamp2(c(-3, 0, 3), c("blue","white", "red")), show_row_names = FALSE, cluster_rows = FALSE)

# Draw the combined heatmap
png(paste(heatmap_name, "colorful", ".png", sep = ""), height = 2000, width = 2000, res = 300)
draw(heatmap_log2_comp_16_vs_rt_2 +
    heatmap_hclust_rep_rt_1 +
    heatmap_hclust_rep_rt_2 +
    heatmap_hclust_rep_rt_3 +
    heatmap_hclust_rep_rt_4 +
    heatmap_hclust_rep_rt_5 +
    heatmap_hclust_rep_16_1 +
    heatmap_hclust_rep_16_2 +
    heatmap_hclust_rep_16_3 +
    heatmap_hclust_rep_16_4 +
    heatmap_hclust_rep_16_5
)
dev.off()

# Save the heatmaps as a pdf
pdf(paste(heatmap_name, "colorful", ".pdf", sep = ""), height = 8, width = 8)
draw(heatmap_log2_comp_16_vs_rt_2 +
    heatmap_hclust_rep_rt_1 +
    heatmap_hclust_rep_rt_2 +
    heatmap_hclust_rep_rt_3 +
    heatmap_hclust_rep_rt_4 +
    heatmap_hclust_rep_rt_5 +
    heatmap_hclust_rep_16_1 +
    heatmap_hclust_rep_16_2 +
    heatmap_hclust_rep_16_3 +
    heatmap_hclust_rep_16_4 +
    heatmap_hclust_rep_16_5
)
dev.off()