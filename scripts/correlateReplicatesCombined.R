# Load the necessary libraries
library(ggplot2)
library(reshape2)

# Input data
args <- commandArgs(trailingOnly = TRUE)
level <- as.character(args[3])

# Read the matrix from a file
mat1 <- as.matrix(read.table(args[1], header = TRUE, check.names = FALSE))
mat2 <- as.matrix(read.table(args[2], header = TRUE, check.names = FALSE))
mat <- cbind(mat1,mat2)
head(mat)

# Take the log of each element in the matrix, adding 1 to each element first
dataset <- log(mat + 1)

cormat <- round(cor(dataset), 3)
melted_cormat <- melt(cormat)

# Make heatmap
# png(paste("correlation_heatmap_", level, "_combined.png", sep = ""), height = 2000, width = 2000, res = 300)

heatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "yellow",
    high = "blue",
    mid = "white",
    midpoint = 0.75,
    limit = c(0.5, 1),
    space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_fixed()

# Save the heatmap as png and PDF files
ggsave(paste("correlation_heatmap_", level, "_combined.png", sep = ""), heatmap, width = 2000, height = 2000, units = "px", dpi = 300)
ggsave(paste("correlation_heatmap_", level, "_combined.pdf", sep = ""), heatmap, width = 2000, height = 2000, units = "px", dpi = 300)