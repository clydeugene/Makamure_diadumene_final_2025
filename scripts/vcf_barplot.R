# Load libraries
library(tidyverse)
library(broom)

# Load data
cliprmvd_readcount <- as.data.frame(read_tsv("data/processed/gatk/genome/cliprmvd_readcounts.tsv", col_names = FALSE))

downsampled_readcount <- as.data.frame(read_tsv("data/processed/gatk/genome/downsampled_readcounts.tsv", col_names = FALSE))

# Add columns names
colnames(cliprmvd_readcount) <- c("sample", "read_count")
colnames(downsampled_readcount)  <- c("sample", "read_count")

# Add a column to indicate whether the sample was downsampled or not
cliprmvd_readcount$downsampled <- "No"
downsampled_readcount$downsampled <- "Yes"


# Combine the cliprmvd and downsampled read counts such that each sample has two rows and a column indicating whether the sample was downsampled or not
combined_readcount <- rbind(cliprmvd_readcount, downsampled_readcount)

# Add a grouping column that indicates whether the sample is from the room temperature or low temperature group
combined_readcount <- combined_readcount %>%
  mutate(temperature = ifelse(grepl("rt", sample), "Room", "Low"))