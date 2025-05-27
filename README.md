# VCF Analysis Pipeline Documentation

## Table of Contents
1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Output Description](#output-description)
6. [Troubleshooting](#troubleshooting)

## Overview

This pipeline performs comprehensive analysis of VCF (Variant Call Format) files, generating detailed reports and statistics for variant calling quality assessment. It's particularly designed for analyzing Plasmodium falciparum whole genome sequencing data.

## Prerequisites

The following software and tools are required:

1. **BCFtools** (version 1.9 or higher)
   - Required for VCF file manipulation and basic statistics
   - [Download and Installation Guide](https://samtools.github.io/bcftools/)

2. **VCFtools** (version 0.1.16 or higher)
   - Required for advanced VCF analysis
   - [Download and Installation Guide](https://vcftools.github.io/downloads.html)

3. **R** (version 4.0.0 or higher) with the following packages:
   - ggplot2
   - dplyr
   - tidyr
   - knitr
   - viridis
   - gridExtra

## Installation

### Option 1: Manual Installation

1. **Install BCFtools:**
   ```bash
   # For Ubuntu/Debian
   sudo apt-get install bcftools

   # For CentOS/RHEL
   sudo yum install bcftools

   # For macOS
   brew install bcftools
   ```

2. **Install VCFtools:**
   ```bash
   # For Ubuntu/Debian
   sudo apt-get install vcftools

   # For CentOS/RHEL
   sudo yum install vcftools

   # For macOS
   brew install vcftools
   ```

3. **Install R and Required Packages:**
   ```bash
   # Install R packages
   R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'knitr', 'viridis', 'gridExtra'), repos='https://cran.rstudio.com/')"
   ```

### Option 2: Using Conda (Recommended)

1. **Install Miniconda** if not already installed
2. **Create and activate environment:**
   ```bash
   conda create -n vcf_analysis -c bioconda -c conda-forge bcftools vcftools r-essentials r-ggplot2 r-dplyr r-tidyr r-knitr r-viridis r-gridextra
   conda activate vcf_analysis
   ```

## Usage

1. **Prepare Your Data:**
   - Ensure your VCF file is compressed (`.vcf.gz`) and indexed (`.vcf.gz.tbi`)
   - Place your VCF file in the working directory

2. **Run the Analysis:**
   ```bash
   # Make the script executable
   chmod +x run_vcf_analysis.sh

   # Run the analysis
   ./run_vcf_analysis.sh
   ```

3. **View Results:**
   - All outputs will be stored in the `VCF_Analysis_Outputs` directory
   - The main report will be available as `vcf_analysis_report.html`

## Output Description

The pipeline generates several output files in the `VCF_Analysis_Outputs` directory:

1. **Basic Statistics:**
   - `snps_only.vcf.gz`: Filtered VCF containing only SNPs
   - `genotype_table.csv`: Genotype information for all variants
   - `allelic_depth.csv`: Allele depth information
   - `depth_per_site.csv`: Read depth per site

2. **Quality Metrics:**
   - `site_qualities.csv`: Quality scores for each variant
   - `allele_frequencies.frq`: Allele frequency statistics
   - `mean_depth_individual.idepth`: Mean depth per individual
   - `mean_depth_site.ldepth.mean`: Mean depth per site

3. **Missing Data Analysis:**
   - `missing_individuals.imiss`: Missing data per individual
   - `missing_sites.lmiss`: Missing data per site

4. **Heterozygosity Analysis:**
   - `heterozygosity.het`: Heterozygosity and inbreeding coefficients

5. **Report:**
   - `vcf_analysis_report.html`: Comprehensive HTML report with visualizations

## Troubleshooting

### Common Issues and Solutions:

1. **BCFtools/VCFtools not found:**
   ```bash
   # Check if tools are in PATH
   which bcftools
   which vcftools
   
   # If not found, add to PATH or use full path
   export PATH=/path/to/tools:$PATH
   ```

2. **R package installation fails:**
   ```bash
   # Try installing packages individually
   R -e "install.packages('ggplot2', repos='https://cran.rstudio.com/')"
   ```

3. **Permission denied:**
   ```bash
   # Make script executable
   chmod +x run_vcf_analysis.sh
   ```

4. **Memory issues:**
   - Ensure sufficient RAM (recommended: 16GB minimum)
   - For large VCF files, consider increasing system swap space

### Getting Help:
- For BCFtools issues: [BCFtools GitHub](https://github.com/samtools/bcftools)
- For VCFtools issues: [VCFtools GitHub](https://github.com/vcftools/vcftools)
- For R package issues: [R Documentation](https://www.r-project.org/other-docs.html)

# VCF Analysis Report Documentation

This documentation explains each plot and statistic in the VCF Analysis Report, providing guidance on interpretation and the significance of each metric in the context of variant calling analysis.

## 1. Sample Information

**Purpose:** This section lists all samples included in the VCF file.

**Interpretation:** The table displays sample IDs, which are important for tracking which individuals are included in the analysis. The number of samples influences statistical power and representation.

**Importance:** Sample information helps verify that all expected samples are present and correctly labeled, which is crucial for downstream analyses, especially in population genetics studies.

## 2. Variant Quality Metrics

### Distribution of Variant Quality Scores

**Purpose:** This histogram shows the distribution of QUAL scores across all variants in the VCF file.

**Interpretation:** 
- The x-axis (log scale) shows quality scores, while the y-axis shows the count of variants.
- Higher quality scores indicate greater confidence in variant calls.
- Ideally, you want to see most variants having high quality scores.
- A bimodal distribution may indicate two distinct classes of variants (high confidence vs. low confidence).

**Importance:** Quality scores help identify potential false positives. Variants with very low quality scores may be filtered out in downstream analyses to improve reliability.

### Quality Score Distribution (Boxplot)

**Purpose:** This boxplot summarizes the central tendency and spread of quality scores.

**Interpretation:**
- The box represents the interquartile range (IQR, 25th to 75th percentile).
- The line inside the box is the median.
- Whiskers typically extend to 1.5 times the IQR.
- Outliers are plotted as individual points.

**Importance:** The boxplot provides a quick summary of quality score distribution, helping identify if most variants are of acceptable quality and detecting outliers.

## 3. Depth Analysis

### Mean Depth per Sample

**Purpose:** This bar chart shows the average sequencing depth (coverage) for each sample.

**Interpretation:**
- The x-axis shows sample IDs, while the y-axis shows mean depth.
- Higher values indicate more sequencing data available for that sample.
- Ideally, all samples should have similar depths to ensure comparable variant calling sensitivity.
- Samples with very low depth may have more false negatives due to insufficient coverage.

**Importance:** Sequencing depth directly impacts the ability to accurately call variants. Samples with unusually low depth may need to be excluded or analyzed with caution.

### Distribution of Mean Depth Across Sites

**Purpose:** This density plot shows the distribution of sequencing depth across genomic positions.

**Interpretation:**
- The x-axis shows depth values, while the y-axis shows density.
- The shape indicates how evenly coverage is distributed across the genome.
- Multiple peaks may indicate systematic biases in coverage.
- A long right tail suggests some regions have extremely high coverage.

**Importance:** Regions with extremely low coverage may have poor variant calling sensitivity, while regions with unusually high coverage could indicate mapping artifacts or copy number variations.

## 4. Missing Data Analysis

### Missing Data by Sample

**Purpose:** This bar chart shows the proportion of missing genotypes for each sample.

**Interpretation:**
- The x-axis shows sample IDs, while the y-axis shows the fraction of missing genotypes.
- Lower values are better, indicating more complete data.
- Samples with high missingness (e.g., >10-20%) may need to be excluded.
- Uneven missingness across samples can bias population genetic analyses.

**Importance:** Missing data can impact many analyses, especially those involving population structure or relatedness. Samples with excessive missing data may need to be removed.

### Distribution of Missing Data Across Sites

**Purpose:** This histogram shows how missingness is distributed across genomic sites.

**Interpretation:**
- The x-axis shows the fraction of samples missing data at a given position, while the y-axis shows the count of sites.
- Ideally, most sites should have low missingness.
- Sites with high missingness across many samples may be in difficult-to-sequence regions.

**Importance:** Sites with high missingness across samples may need to be filtered out before certain analyses, as they can reduce statistical power and introduce biases.

## 5. Heterozygosity Analysis

### Inbreeding Coefficient by Sample

**Purpose:** This bar chart shows the inbreeding coefficient (F) for each sample.

**Interpretation:**
- The x-axis shows sample IDs, while the y-axis shows the inbreeding coefficient.
- F values around 0 suggest random mating (Hardy-Weinberg equilibrium).
- Positive F values indicate an excess of homozygosity (inbreeding).
- Negative F values suggest an excess of heterozygosity.
- For haploid organisms like Plasmodium, high heterozygosity may indicate mixed infections.

**Importance:** The inbreeding coefficient can reveal population structure, inbreeding, sample contamination, or mixed infections in parasites.

### Observed vs Expected Homozygosity

**Purpose:** This scatter plot compares observed versus expected homozygous sites for each sample.

**Interpretation:**
- Points along the diagonal indicate samples in Hardy-Weinberg equilibrium.
- Points above the diagonal have more homozygosity than expected (potential inbreeding).
- Points below the diagonal have more heterozygosity than expected (potential outbreeding or sample contamination).

**Importance:** Deviations from expected heterozygosity can indicate sample quality issues, population structure, or selection.

## 6. Detailed Summary Statistics

### Quality Statistics

**Purpose:** This table summarizes key metrics about variant quality scores.

**Interpretation:**
- Total Sites: The total number of variant sites analyzed.
- Mean and Median Quality: Central tendency measures of quality scores.
- Quality SD: Indicates the spread of quality scores.
- Quality Range: Shows the minimum and maximum quality scores.

**Importance:** These statistics provide a quick overview of variant quality, helping identify potential issues with the variant calling process.

### Depth Statistics

**Purpose:** This table summarizes sequencing depth across the dataset.

**Interpretation:**
- Mean Depth: Average coverage across all samples and sites.
- Depth SD: Variation in depth, with high values indicating uneven coverage.
- Min/Max Depth: Extremes of coverage.
- Sites < 10x: Count of sites with potentially insufficient coverage for reliable variant calling.

**Importance:** Depth statistics help evaluate the overall quality of sequencing data and identify potential limitations in variant detection.

### Missing Data Statistics

**Purpose:** This table summarizes missing data across the dataset.

**Interpretation:**
- Mean Missing: Average proportion of missing data across all samples.
- Max Missing: Highest proportion of missing data in any sample.
- Samples > 20% Missing: Count of samples with potentially problematic levels of missing data.
- Total Missing Sites: Sum of all missing genotypes across all samples and sites.

**Importance:** These statistics help evaluate data completeness and identify patterns of missingness that might affect analyses.

## 7. Allele Frequency Analysis

### Allele Frequency Spectrum

**Purpose:** This histogram shows the distribution of allele frequencies across all variable sites.

**Interpretation:**
- The x-axis shows allele frequency, while the y-axis shows the count of variants.
- The shape of this distribution can reveal population history and selection.
- For neutral evolution, many rare variants and fewer common variants are expected.
- Peaks at intermediate frequencies may indicate balancing selection.
- For Plasmodium, the shape can indicate drug resistance patterns or adaptation.

**Importance:** The allele frequency spectrum is fundamental for population genetics analyses and can reveal evolutionary processes and selective pressures affecting the population.

## Conclusion

This VCF analysis report provides a comprehensive quality assessment of variant calls, helping identify potential issues with sequencing data, variant calling, or sample quality. Key areas to examine include:

1. **Sample quality:** Check for samples with unusual depth, missing data, or heterozygosity.
2. **Variant quality:** Examine the distribution of quality scores to identify potential false positives.
3. **Coverage uniformity:** Assess whether depth is consistent across samples and genomic regions.
4. **Completeness:** Evaluate the extent of missing data and its distribution.
5. **Population genetics indicators:** Use heterozygosity and allele frequency metrics to understand the genetic structure of the population.

By addressing potential issues identified in this report, you can improve the reliability of downstream analyses and interpretations of genetic variation in your samples. 
