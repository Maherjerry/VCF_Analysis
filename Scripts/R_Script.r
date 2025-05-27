# VCF Analysis Report Generator
# Author: Bah Muhammed
# Description: Generates analysis plots and statistics from VCF analysis outputs

# Load required libraries
required_packages <- c("ggplot2", "dplyr", "tidyr", "knitr", "viridis", "gridExtra")
for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, repos = "https://cran.rstudio.com/")
    library(package, character.only = TRUE)
  }
}

# Set working directory to the output directory
output_dir <- "VCF_Analysis_Outputs"
if (!dir.exists(output_dir)) {
  stop("Output directory does not exist. Please run the VCF analysis script first.")
}
setwd(output_dir)

# Create plots directory
plots_dir <- "plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Function to save plots
save_plot <- function(plot, filename, width = 12, height = 6) {
  pdf(file.path(plots_dir, filename), width = width, height = height)
  print(plot)
  dev.off()
}

# Read and process data
tryCatch({
  # Sample information
  samples <- read.csv("sample_names.csv", header = FALSE)
  colnames(samples) <- c("Sample")
  
  # Quality metrics
  qualities <- read.csv("site_qualities.csv", sep = "\t", header = FALSE)
  colnames(qualities) <- c("CHROM", "POS", "REF", "ALT", "QUAL")
  qualities$QUAL <- as.numeric(qualities$QUAL)
  
  # Depth analysis
  indv_depth <- read.csv("mean_depth_individual.idepth", sep = "\t")
  mean_site_depth <- read.csv("mean_depth_site.ldepth.mean", sep = "\t")
  
  # Missing data analysis
  missing_indv <- read.csv("missing_individuals.imiss", sep = "\t")
  missing_sites <- read.csv("missing_sites.lmiss", sep = "\t")
  
  # Heterozygosity analysis
  het <- read.csv("heterozygosity.het", sep = "\t")
  
  # Generate plots
  
  # 1. Quality Score Distribution
  p1 <- ggplot(qualities, aes(x = QUAL)) + 
    geom_histogram(bins = 50, fill = "steelblue") +
    theme_minimal() +
    labs(title = "Distribution of Variant Quality Scores",
         x = "Quality Score", y = "Count") +
    scale_x_continuous(trans = "log10")
  
  p2 <- ggplot(qualities, aes(y = QUAL)) +
    geom_boxplot(fill = "lightblue") +
    theme_minimal() +
    labs(title = "Quality Score Distribution",
         y = "Quality Score")
  
  save_plot(gridExtra::grid.arrange(p1, p2, ncol = 2), "quality_metrics.pdf")
  
  # 2. Depth Analysis
  p1 <- ggplot(indv_depth, aes(x = INDV, y = MEAN_DEPTH)) +
    geom_col(fill = "forestgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Mean Depth per Sample", 
         x = "Sample", y = "Mean Depth")
  
  p2 <- ggplot(mean_site_depth, aes(x = MEAN_DEPTH)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    theme_minimal() +
    labs(title = "Distribution of Mean Depth Across Sites",
         x = "Mean Depth", y = "Density")
  
  save_plot(gridExtra::grid.arrange(p1, p2, ncol = 2), "depth_analysis.pdf")
  
  # 3. Missing Data Analysis
  p1 <- ggplot(missing_indv, aes(x = INDV, y = F_MISS)) +
    geom_col(fill = "coral") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Missing Data by Sample", 
         x = "Sample", y = "Fraction of Missing Genotypes")
  
  p2 <- ggplot(missing_sites, aes(x = F_MISS)) +
    geom_histogram(fill = "coral", alpha = 0.5, bins = 50) +
    theme_minimal() +
    labs(title = "Distribution of Missing Data Across Sites",
         x = "Fraction of Missing Genotypes", y = "Count")
  
  save_plot(gridExtra::grid.arrange(p1, p2, ncol = 2), "missing_data_analysis.pdf")
  
  # 4. Heterozygosity Analysis
  p1 <- ggplot(het, aes(x = INDV, y = F)) +
    geom_col(fill = "purple") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Inbreeding Coefficient by Sample", 
         x = "Sample", y = "Inbreeding Coefficient")
  
  p2 <- ggplot(het, aes(x = O.HOM., y = E.HOM.)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = "Observed vs Expected Homozygosity",
         x = "Observed Homozygous Sites", y = "Expected Homozygous Sites")
  
  save_plot(gridExtra::grid.arrange(p1, p2, ncol = 2), "heterozygosity_analysis.pdf")
  
  # Generate summary statistics
  summary_stats <- list(
    Quality = data.frame(
      Metric = c("Total Sites", "Mean Quality", "Median Quality", "Quality SD", "Quality Range"),
      Value = c(
        nrow(qualities),
        round(mean(qualities$QUAL, na.rm = TRUE), 2),
        round(median(qualities$QUAL, na.rm = TRUE), 2),
        round(sd(qualities$QUAL, na.rm = TRUE), 2),
        paste(round(min(qualities$QUAL, na.rm = TRUE), 2), "-", 
              round(max(qualities$QUAL, na.rm = TRUE), 2))
      )
    ),
    Depth = data.frame(
      Metric = c("Mean Depth", "Depth SD", "Min Depth", "Max Depth", "Sites < 10x"),
      Value = c(
        round(mean(indv_depth$MEAN_DEPTH), 2),
        round(sd(indv_depth$MEAN_DEPTH), 2),
        round(min(indv_depth$MEAN_DEPTH), 2),
        round(max(indv_depth$MEAN_DEPTH), 2),
        sum(mean_site_depth$MEAN_DEPTH < 10)
      )
    ),
    Missing = data.frame(
      Metric = c("Mean Missing", "Max Missing", "Samples > 20% Missing", "Total Missing Sites"),
      Value = c(
        round(mean(missing_indv$F_MISS), 4),
        round(max(missing_indv$F_MISS), 4),
        sum(missing_indv$F_MISS > 0.2),
        sum(missing_sites$N_MISS)
      )
    )
  )
  
  # Save summary statistics
  for (stat_name in names(summary_stats)) {
    write.csv(summary_stats[[stat_name]], 
              file.path(plots_dir, paste0(tolower(stat_name), "_statistics.csv")), 
              row.names = FALSE)
  }
  
  cat("Analysis complete. Results saved in", output_dir, "directory.\n")
  cat("Plots saved in", file.path(output_dir, plots_dir), "directory.\n")
  
}, error = function(e) {
  cat("Error occurred:", conditionMessage(e), "\n")
  quit(status = 1)
})