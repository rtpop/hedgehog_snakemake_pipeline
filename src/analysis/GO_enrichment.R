## loading libraries
required_libraries <- c("optparse",
                        "data.table",
                        "topGO")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##

option_list <- list(
    optparse::make_option(c("-g", "--gmt_file"), type = "character", help = "GMT file with the communities."),
    optparse::make_option(c("-b", "--bg_file"), type = "character", help = "Path to a file with the background genes. Only needed if auto_bg is FALSE."),
    optparse::make_option(c("-a", "--auto_bg"), type = "logical", default = TRUE, help = "If TRUE, the background is automatically set to the genes in the gene sets.")
    optparse::make_option(c("-s", "--save_all"), type = "logical", default = TRUE, help = "If TRUE, the GO enrichment for each community will be saved in addition to a summary for all communities.")
    optparse::make_option(c("-t", "--sig_thresh"), type = "numeric", default = 0.05, help = "P-value significant threshold")
    optparse::make_option(c("-st", "--statistic"), type = "character", default = "fisher", help = "Test statistic to be used.")
    optparse::make_option(c("-ag", "--algorithm"), type = "character", default = "classic", help = "topGO algorithm to use for enrichment.")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## --------- ##
## Variables ##
## --------- ##

GMT_FILE <- opt$gmt_file
BG_FILE <- opt$bg_file
AUTO_BG <- opt$auto_bg
SAVE <- opt$save_all
THRESHOLD <- opt$sig_thresh
STAT <- opt$statistic
ALG <- opt$algorithm

## source functions
source("src/analysis/GO_enrichment_fn.R")

## Run
run_go_enrichment(gmt_file = GMT_FILE, bg_file = BG_FILE, save_all = SAVE, sig_thresh = THRESHOLD, statistic = STAT, algorithm = ALG)