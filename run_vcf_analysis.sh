#!/bin/bash

# Set input and output paths
VCF_FILE="1332-PF-SL-NGWA-SM_genotyped.vcf.gz"
OUTPUT_DIR="VCF_Analysis_Outputs"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "Starting VCF summary analysis..."

##########################
# BCFTOOLS-based exports
##########################

# Function to check if a file exists and is not empty
check_file() {
    if [ -f "$1" ] && [ -s "$1" ]; then
        return 0  # File exists and is not empty
    else
        return 1  # File doesn't exist or is empty
    fi
}

# Function to run a command if output doesn't exist
run_if_missing() {
    local output_file=$1
    local command=$2
    local description=$3

    if ! check_file "$output_file"; then
        echo "$description..."
        eval "$command"
    else
        echo "Skipping $description (output already exists)"
    fi
}

# BCFTOOLS commands
run_if_missing "${OUTPUT_DIR}/snps_only.vcf.gz" \
    "bcftools view -v snps -Oz -o ${OUTPUT_DIR}/snps_only.vcf.gz ${VCF_FILE} && bcftools index ${OUTPUT_DIR}/snps_only.vcf.gz" \
    "Extracting SNPs"

run_if_missing "${OUTPUT_DIR}/genotype_table.csv" \
    "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${VCF_FILE} > ${OUTPUT_DIR}/genotype_table.csv" \
    "Extracting genotype table"

run_if_missing "${OUTPUT_DIR}/allelic_depth.csv" \
    "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${VCF_FILE} > ${OUTPUT_DIR}/allelic_depth.csv" \
    "Extracting allelic depth"

run_if_missing "${OUTPUT_DIR}/depth_per_site.csv" \
    "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' ${VCF_FILE} > ${OUTPUT_DIR}/depth_per_site.csv" \
    "Extracting read depth"

run_if_missing "${OUTPUT_DIR}/site_qualities.csv" \
    "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' ${VCF_FILE} > ${OUTPUT_DIR}/site_qualities.csv" \
    "Extracting site quality"

run_if_missing "${OUTPUT_DIR}/sample_names.csv" \
    "bcftools query -l ${VCF_FILE} > ${OUTPUT_DIR}/sample_names.csv" \
    "Extracting sample names"

##########################
# VCFTOOLS-based exports
##########################

run_if_missing "${OUTPUT_DIR}/allele_frequencies.frq" \
    "vcftools --gzvcf ${VCF_FILE} --freq2 --out ${OUTPUT_DIR}/allele_frequencies --max-alleles 2" \
    "Calculating allele frequency"

run_if_missing "${OUTPUT_DIR}/mean_depth_individual.idepth" \
    "vcftools --gzvcf ${VCF_FILE} --depth --out ${OUTPUT_DIR}/mean_depth_individual" \
    "Calculating mean depth per individual"

run_if_missing "${OUTPUT_DIR}/mean_depth_site.ldepth.mean" \
    "vcftools --gzvcf ${VCF_FILE} --site-mean-depth --out ${OUTPUT_DIR}/mean_depth_site" \
    "Calculating mean depth per site"

run_if_missing "${OUTPUT_DIR}/missing_individuals.imiss" \
    "vcftools --gzvcf ${VCF_FILE} --missing-indv --out ${OUTPUT_DIR}/missing_individuals" \
    "Calculating missing data per individual"

run_if_missing "${OUTPUT_DIR}/missing_sites.lmiss" \
    "vcftools --gzvcf ${VCF_FILE} --missing-site --out ${OUTPUT_DIR}/missing_sites" \
    "Calculating missing data per site"

run_if_missing "${OUTPUT_DIR}/heterozygosity.het" \
    "vcftools --gzvcf ${VCF_FILE} --het --out ${OUTPUT_DIR}/heterozygosity" \
    "Calculating heterozygosity and inbreeding coefficient"

##########################
# Generate R Markdown Report
##########################

# Check if we need to regenerate the report
if [ ! -f "${OUTPUT_DIR}/vcf_analysis_report.html" ] || [ "$(find ${OUTPUT_DIR} -name "*.csv" -o -name "*.imiss" -o -name "*.lmiss" -o -name "*.het" -o -name "*.frq" -o -name "*.idepth" -o -name "*.ldepth.mean" -newer "${OUTPUT_DIR}/vcf_analysis_report.html" 2>/dev/null)" ]; then
    echo "Generating R Markdown report..."
    
    # Create R Markdown report
    cat > ${OUTPUT_DIR}/vcf_report.Rmd << 'EOL'
---
title: "VCF Analysis Report"
author: "Plasmodium falciparum WGS Analysis Pipeline"
date: "`r Sys.Date()`"
output: 
  html_document:
    theme: cosmo
    toc: true
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(viridis)
library(gridExtra)
```

## Overview

This report summarizes the variant calls from Plasmodium falciparum whole genome sequencing.

## Sample Information

```{r}
samples <- read.csv("sample_names.csv", header=FALSE)
colnames(samples) <- c("Sample")
kable(samples, caption = "List of Samples")
```

## Variant Quality Metrics

```{r}
qualities <- read.csv("site_qualities.csv", sep="\t", header=FALSE)
colnames(qualities) <- c("CHROM", "POS", "REF", "ALT", "QUAL")
qualities$QUAL <- as.numeric(qualities$QUAL)

# Quality Score Distribution
p1 <- ggplot(qualities, aes(x=QUAL)) + 
  geom_histogram(bins=50, fill="steelblue") +
  theme_minimal() +
  labs(title="Distribution of Variant Quality Scores",
       x="Quality Score", y="Count") +
  scale_x_continuous(trans="log10")

# Quality Score Summary
p2 <- ggplot(qualities, aes(y=QUAL)) +
  geom_boxplot(fill="lightblue") +
  theme_minimal() +
  labs(title="Quality Score Distribution",
       y="Quality Score")

# Display plots side by side
gridExtra::grid.arrange(p1, p2, ncol=2)
```

## Depth Analysis

```{r}
indv_depth <- read.csv("mean_depth_individual.idepth", sep="\t")
mean_site_depth <- read.csv("mean_depth_site.ldepth.mean", sep="\t")

# Mean Depth per Sample
p1 <- ggplot(indv_depth, aes(x=INDV, y=MEAN_DEPTH)) +
  geom_col(fill="forestgreen") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Mean Depth per Sample", 
       x="Sample", y="Mean Depth")

# Site Depth Distribution
p2 <- ggplot(mean_site_depth, aes(x=MEAN_DEPTH)) +
  geom_density(fill="lightblue", alpha=0.5) +
  theme_minimal() +
  labs(title="Distribution of Mean Depth Across Sites",
       x="Mean Depth", y="Density")

# Display plots side by side
gridExtra::grid.arrange(p1, p2, ncol=2)
```

## Missing Data Analysis

```{r}
missing_indv <- read.csv("missing_individuals.imiss", sep="\t")
missing_sites <- read.csv("missing_sites.lmiss", sep="\t")

# Missing Data by Sample
p1 <- ggplot(missing_indv, aes(x=INDV, y=F_MISS)) +
  geom_col(fill="coral") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Missing Data by Sample", 
       x="Sample", y="Fraction of Missing Genotypes")

# Missing Data Distribution
p2 <- ggplot(missing_sites, aes(x=F_MISS)) +
  geom_histogram(fill="coral", alpha=0.5, bins=50) +
  theme_minimal() +
  labs(title="Distribution of Missing Data Across Sites",
       x="Fraction of Missing Genotypes", y="Count")

# Display plots side by side
gridExtra::grid.arrange(p1, p2, ncol=2)
```

## Heterozygosity Analysis

```{r}
het <- read.csv("heterozygosity.het", sep="\t")

# Inbreeding Coefficient
p1 <- ggplot(het, aes(x=INDV, y=F)) +
  geom_col(fill="purple") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Inbreeding Coefficient by Sample", 
       x="Sample", y="Inbreeding Coefficient")

# Heterozygosity Distribution
p2 <- ggplot(het, aes(x=O.HOM., y=E.HOM.)) +
  geom_point(alpha=0.5) +
  theme_minimal() +
  labs(title="Observed vs Expected Homozygosity",
       x="Observed Homozygous Sites", y="Expected Homozygous Sites")

# Display plots side by side
gridExtra::grid.arrange(p1, p2, ncol=2)
```

## Detailed Summary Statistics

```{r}
# Quality Statistics
quality_stats <- data.frame(
  Metric = c("Total Sites", "Mean Quality", "Median Quality", "Quality SD", "Quality Range"),
  Value = c(
    nrow(qualities),
    round(mean(qualities$QUAL, na.rm=TRUE), 2),
    round(median(qualities$QUAL, na.rm=TRUE), 2),
    round(sd(qualities$QUAL, na.rm=TRUE), 2),
    paste(round(min(qualities$QUAL, na.rm=TRUE), 2), "-", round(max(qualities$QUAL, na.rm=TRUE), 2))
  )
)

# Depth Statistics
depth_stats <- data.frame(
  Metric = c("Mean Depth", "Depth SD", "Min Depth", "Max Depth", "Sites < 10x"),
  Value = c(
    round(mean(indv_depth$MEAN_DEPTH), 2),
    round(sd(indv_depth$MEAN_DEPTH), 2),
    round(min(indv_depth$MEAN_DEPTH), 2),
    round(max(indv_depth$MEAN_DEPTH), 2),
    sum(mean_site_depth$MEAN_DEPTH < 10)
  )
)

# Missing Data Statistics
missing_stats <- data.frame(
  Metric = c("Mean Missing", "Max Missing", "Samples > 20% Missing", "Total Missing Sites"),
  Value = c(
    round(mean(missing_indv$F_MISS), 4),
    round(max(missing_indv$F_MISS), 4),
    sum(missing_indv$F_MISS > 0.2),
    sum(missing_sites$N_MISS)
  )
)

# Display all statistics
kable(quality_stats, caption = "Quality Statistics")
kable(depth_stats, caption = "Depth Statistics")
kable(missing_stats, caption = "Missing Data Statistics")
```

## Allele Frequency Analysis

```{r}
# Read allele frequencies with row.names=NULL to prevent duplicate row names
allele_freq <- read.csv("allele_frequencies.frq", sep="\t", check.names=FALSE, row.names=NULL)

# Print column names for debugging
print("Column names in allele frequency file:")
print(colnames(allele_freq))

# Allele Frequency Spectrum
ggplot(allele_freq, aes(x=`{FREQ}`)) +
  geom_histogram(bins=50, fill="purple", alpha=0.5) +
  theme_minimal() +
  labs(title="Allele Frequency Spectrum",
       x="Allele Frequency", y="Count")
```
EOL

    # Change to output directory for R rendering
    cd ${OUTPUT_DIR}

    # Install required R packages if not already installed
    Rscript -e "if(!require('gridExtra')) install.packages('gridExtra', repos='http://cran.rstudio.com/')"
    Rscript -e "if(!require('viridis')) install.packages('viridis', repos='http://cran.rstudio.com/')"

    # Render the report
    Rscript -e "rmarkdown::render('vcf_report.Rmd', output_file = 'vcf_analysis_report.html')"
else
    echo "Skipping report generation (report is up to date)"
fi

echo "Analysis complete! Results are in ${OUTPUT_DIR}"