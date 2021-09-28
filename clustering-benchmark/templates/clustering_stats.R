#!/usr/bin/env Rscript --vanilla

# Script to 
suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(readr)
    library(dplyr)
})

# load all clustering summaries
result_files <- list.files(".", pattern = "*_psm.tsv\$", full.names = TRUE)

if (length(result_files) < 1) {
    stop("Error: Failed to find any result files.")
}

message("Found ", length(result_files), " clustering results")

raw_data <- do.call(bind_rows, lapply(result_files, function(filename) {
    suppressMessages({
        data <- read_tsv(filename)
    })

    data\$tool <- gsub(".*/(.*)_psm.tsv", "\\\\1", filename)

    return(data)
}))

raw_data <- raw_data %>%
    mutate(cluster_id = paste0(tool, "_", cluster))

# get the most common peptide petide per cluster
max_peptide <- raw_data %>%
    filter(!is.na(peptide)) %>%
    group_by(cluster_id, peptide) %>%
    summarise(n_psms_max_peptide = n()) %>%
    arrange(desc(n_psms_max_peptide)) %>%
    filter(!duplicated(cluster_id)) %>%
    rename(max_peptide = peptide)

# merge on the cluster level
cluster_data <- raw_data %>%
    left_join(max_peptide, by = "cluster_id") %>%
    group_by(tool, cluster) %>%
    summarise(
        n_spectra = n(),
        n_psms = sum(!is.na(peptide)),
        max_peptide = unique(max_peptide),
        max_ratio = unique(n_psms_max_peptide) / sum(!is.na(peptide))
    ) %>%
    mutate(max_ratio = ifelse(n_psms > 0, max_ratio, NA))

# get the stats per tool
tool_stats <- cluster_data %>%
    filter(n_spectra >= 3) %>%
    group_by(tool) %>%
    summarise(
        clustered_spectra = sum(n_spectra), 
        clustered_psms = sum(n_psms), 
        rel_incorrect_spectra = sum(n_psms * (1 - max_ratio), na.rm = T) / sum(n_psms, na.rm = T))

write_tsv(tool_stats, file = "tool_stats.tsv")
 
# create the summary plot
n_spectra <- nrow(raw_data %>%
    filter(tool == tool[1]))

plot_data <- tool_stats %>%
    mutate(
        tool_name = gsub("_.*", "", tool), 
        tool_settings = gsub(".*_", "", tool),
        rel_clustered_spectra = clustered_spectra / n_spectra)

summary_plot <- ggplot(plot_data, aes(x = rel_incorrect_spectra, y = rel_clustered_spectra, group = tool_name, color = tool_name)) +
    geom_point() +
    geom_line() +
    scale_color_brewer(palette = "Set1")

png(width = 1500, height = 1000, res = 150)
print(summary_plot)
dev.off()
