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
    optparse::make_option(c("-g", "--gmt-file"), type = "character", help = "GMT file with the communities."),
    optparse::make_option(c("-b", "--bg-file"), type = "character", help = "Path to a file with the background genes. Only needed if auto_bg is FALSE."),
    optparse::make_option(c("-a", "--auto-bg"), type = "logical", default = TRUE, help = "If TRUE, the background is automatically set to the genes in the gene sets."),
    optparse::make_option(c("-s", "--save-all"), type = "logical", default = TRUE, help = "If TRUE, the GO enrichment for each community will be saved in addition to a summary for all communities."),
    optparse::make_option(c("-t", "--sig-thresh"), type = "numeric", default = 0.05, help = "P-value significant threshold"),
    optparse::make_option(c("--statistic"), type = "character", default = "fisher", help = "Test statistic to be used."),
    optparse::make_option(c("--algorithm"), type = "character", default = "classic", help = "topGO algorithm to use for enrichment."),
    optparse::make_option(c("-o", "--output-dir"), type = "character", help = "Path to output directory")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## --------- ##
## Variables ##
## --------- ##


GMT_FILE <- opt$`gmt-file`
BG_FILE <- opt$`bg-file`
AUTO_BG <- opt$`auto-bg`
SAVE <- opt$`save-all`
THRESHOLD <- opt$`sig-thresh`
STAT <- opt$statistic
ALG <- opt$algorithm
OUT_DIR <- opt$`output-dir`

## source functions
source("src/analysis/GO_enrichment_fn.R")
source("src/utils/utils.R")

## Run
run_go_enrichment(gmt_file = GMT_FILE, out_dir = OUT_DIR, bg_file = BG_FILE, save_all = SAVE, sig_thresh = THRESHOLD, statistic = STAT, algorithm = ALG)