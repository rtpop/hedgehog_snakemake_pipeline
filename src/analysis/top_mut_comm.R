# load libraries
required_libraries <- c("optparse",
                        "data.table",
                        "dplyr")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##
option_list <- list(
    optparse::make_option(c("-c", "--communities"), type = "character", help = "File path to gmt file with communities."),
    optparse::make_option(c("-s", "--scores"), type = "character", help = "File path for community scores."),
    optparse::make_option(c("-m", "--mut"), type = "character", help = "File path to mutation scores."),
    optparse::make_option(c("-n", "--n-comm"), type = "integer", default = 10, help = "Number of 'top' communities."),
    optparse::make_option(c("-o", "--output-dir"), type = "character", help = "Path to output directory")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
print(opt)

# Variables
GMT_FILE <- opt$communities
COMM_SCORES <- opt$scores
MUT <- opt$mut
N_COMM <- opt$`n-comm`
OUT_DIR <- opt$`output-dir`

# source functions
source("src/analysis/top_mut_comm_fn.R")
source("src/utils/utils.R")

## ---------------- ##
## Select top genes ##
## ---------------- ##
top_mut_comm <- select_top_mut_comm(communities = GMT_FILE, scores = COMM_SCORES, mut = MUT, n_comm = N_COMM)
gene_count_summary <- mut_gene_summary(top_mut_comm)

## Save
save(top_mut_comm, file = file.path(OUT_DIR, paste0("top_mut_comm_", N_COMM, ".RData")))
fwrite(gene_count_summary, file.path(OUT_DIR, paste0("gene_mutation_summary_", N_COMM, ".txt")), sep = "\t")