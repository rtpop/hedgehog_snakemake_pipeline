required_libraries <- c("data.table")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##

option_list <- list(
    optparse::make_option(c("-f", "--files"), type = "character", help = "Comma-separated list of paths to the filtering benchmark result files."),
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to output file")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# parse the input files
input_files <- strsplit(opt$files, ",")[[1]]

## --------- ##
## Variables ##
## --------- ##
OUT <- opt$output

consolidate_data <- function(files) {
    data_list <- lapply(files, function(file) {
        if (file.exists(file)) {
            data <- fread(file, header = TRUE, sep = "\t", data.table = FALSE)
            # Ensure the 'file' column is present
            if (!"file" %in% colnames(data)) {
                data$file <- file
            }
            return(data)
        } else {
            warning(paste("File not found:", file))
            return(NULL)
        }
    })

    # Remove NULL entries (files that were not found)
    data_list <- Filter(Negate(is.null), data_list)

    # Combine all data frames into one
    consolidated_data <- do.call(rbind, data_list)

    # Ensure the 'file' column is a factor
    consolidated_data$file <- as.factor(consolidated_data$file)

    fwrite(consolidated_data, file = OUT, sep = "\t", row.names = FALSE, quote = FALSE)
    message(paste("Consolidated data written to:", OUT))
}

consolidate_data(input_files)
# End of the script