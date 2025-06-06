##----------------------------------------------------------##
## List of R packages requires in the 'minimal_output' mode ##
##----------------------------------------------------------##

## CRAN packages:
required_packages_cran = c(
    "BiocManager",    # To install Bioconductor packages
    # "curl",           # For downloading files
    "dplyr",          # For data manipulation
    "data.table",     # For data manipulation with data.table
    "optparse",       # For reading command-line arguments
    "devtools",       # For installing R packages from GitHub
    "purrr",          # For functional programming tools
    # "rcartocolor",    # For color palettes for maps and plots
    "reshape2",       # For reshaping data frames
    # "this.path",      # For creating relative paths
    "tidyr",          # For tidying data
    "ggplot2",        # For creating plots
    "ggrepel",        # For adding labels to ggplot2 plots
    "ggpubr",         # For ggplot2 publication-ready plots
    "sessioninfo",    # For session information
    "pheatmap",       # For drawing heatmaps
    "RColorBrewer",   # For color palettes for heatmaps
    "stringr",        # For string manipulation
    "magrittr"       # For pipe operator
    # stylo,            # For cosine distance calculation
    # vegan           # For binomial distance calculation
)

install.packages(
    required_packages_cran,
    dependencies = TRUE, verbose = TRUE)

## Bioconductor packages:
required_packages_bioconductor <- c(
    "topGO",          # For Gene Ontology analysis
    "org.Hs.eg.db"    # for GO analysis
)

BiocManager::install(
    required_packages_bioconductor)