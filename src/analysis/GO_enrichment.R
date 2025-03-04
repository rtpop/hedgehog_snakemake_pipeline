## loading libraries
required_libraries <- c("optparse",
                        "data.table",
                        "topGO")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## Command line options
option_list <- list(
    optparse::make_option(c("-g", "--gmt_file"), type = "character", help = "GMT file with the communities."),
    optparse::make_option(c("-b", "--bg_file"), type = "character", help = "Path to a file with the background genes. Only needed if auto_bg is FALSE."),
    optparse::make_option(c("-a", "--auto_bg"), type = "logical", default = TRUE, help = "If TRUE, the background is automatically set to the genes in the gene sets.")
)