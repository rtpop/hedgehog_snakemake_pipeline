## loading libraries
required_libraries <- c("optparse",
                        "data.table",
                        "ggplot2")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##

option_list <- list(
    optparse::make_option(c("-i", "--input"), type = "character", help = "GMT file with the communities."),
    optparse::make_option(c("-s", "--stats"), type = "character", help = "Path to the stats file."),
    optparse::make_option(c("-o", "--out-file"), type = "character", help = "Path to output file"),
    optparse::make_option(c("-p", "--separate"), type = "logical", help = "If plots should be saved separately.")
    
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## --------- ##
## Variables ##
## --------- ##

GMT_FILE <- opt$input
STATS_FILE <- opt$`stats`
OUT_FILE <- opt$`out-file`
TISSUE <- strsplit(STATS_FILE, "/")[[1]][2]
SEPARATE <- opt$`separate`

## source functions
source("src/analysis/community_plot_fn.R")
source("src/utils/utils.R")

## Run
comm <- process_community(gmt_file = GMT_FILE, stats_file = STATS_FILE, tissue = TISSUE)


if (SEPARATE) {
    # the file names for this use case make no sense and will need to edit later
    # Plot number of communities
    plot_n_communities(comm[[1]], TISSUE, file.path(OUT_FILE, "plots", paste0(TISSUE, "_n_communities.png")))
    
    # Plot community size
    plot_community_size(comm[[2]], TISSUE, file.path(OUT_FILE, "plots", paste0(TISSUE, "_community_size.png")))
} else {
    
    # Plot number of communities and community size in one plot
    p <- plot_community_stats(comm[[1]], TISSUE)
    q <- plot_community_size(comm[[2]], TISSUE)

    # Save the combined plot
    ggsave(OUT_FILE, plot = p + q, width = 10, height = 5)

}