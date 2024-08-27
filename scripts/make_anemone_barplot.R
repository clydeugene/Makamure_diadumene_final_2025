# Load the required libraries
library(tidyverse)

# Read in command line arguments. The first five are numeric and the last is a string
args <- commandArgs(trailingOnly = TRUE)
a <- as.numeric(args[1])
b <- as.numeric(args[2])
c <- as.numeric(args[3])
d <- as.numeric(args[4])
e <- as.numeric(args[5])

outfile_name <- as.character(args[6])

# Calculate percentages to display on the barplot
A <- signif(((a + c) / e) * 100, 3)
B <- signif(((b) / e) * 100, 3)
C <- signif((((d - (a + c + b)) + (e - d)) / e) * 100, 3)

# Make a data frame of the percentages
plot_data <- data.frame(
  category = c("Other Anemones and Corals", "Diadumene Genome", "Other IDs and unidentified"),
  percentage = c(A, B, C)
)

# Make a barplot of the percentages
plot_final <- ggplot(plot_data, aes(x = category, y = percentage)) +
  geom_bar(stat = "identity", fill = "gray") +
  geom_text(aes(label = paste0(percentage, "%")), vjust = -0.5) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"))

# Save the barplot as a png, pdf, and svg files
ggsave(paste0(outfile_name, ".png"), plot_final, width = 8, height = 8, dpi = 300)
ggsave(paste0(outfile_name, ".pdf"), plot_final, width = 8, height = 8, dpi = 300)
ggsave(paste0(outfile_name, ".svg"), plot_final, width = 8, height = 8, dpi = 300)


