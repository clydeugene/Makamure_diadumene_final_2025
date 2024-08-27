# Load the necessary libraries
library(ggplot2)
library(reshape2)

# Input data
args <- commandArgs(trailingOnly = TRUE)

# Read the matrix from a file
mat <- as.matrix(read.table(args[1], header = TRUE, check.names = FALSE))

# Read the condition from the command line arguments
condition <- as.character(args[2])



# Take the log of each element in the matrix, adding 1 to each element first
dataset <- log(mat + 1)

cormat <- round(cor(dataset), 3)
melted_cormat <- melt(cormat)

# Make heatmap

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


# Make png and PDF versions of the heatmap
ggsave(paste("correlation_heatmap_", condition, ".png", sep = ""), heatmap, width = 2000, height = 2000, units = "px", dpi = 300)
ggsave(paste("correlation_heatmap_", condition, ".pdf", sep = ""), heatmap, width = 2000, height = 2000, units = "px", dpi = 300)
