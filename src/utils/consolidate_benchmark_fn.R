#' @name prepare_filtering_bench
#' 
#' @title Prepare filtering benchmark results
#' 
#' @description This function takes a data frame containing filtering benchmark results and 
#' adds resolution and tissue type columns to it.
#' 
#' @param file_name File name for data frame containing the filtering benchmark results. The data frame should have the following columns:
#' \itemize{
#'  \item \code{Network}: The method of filtering used (ELAND, Prior or Top).
#'  \item \code{Modularity}: The modularity score of the filtered network.
#'  \item \code{Density}: The density of the filtered network.
#'  \item \code{Number of Edges}: The number of edges in the filtered network.
#' 
#' @param tissue_type The type of tissue used in the filtering benchmark. This should be a string that describes the tissue type.
#' 
#' @return A data frame with the same columns as the input, but with additional columns for resolution and tissue type.
prepare_filtering_bench <- function(file_name, tissue_type) {
    # Read the data from the file
    df <- data.table::fread(file_name, header = TRUE)

    print(file_name)
    
    # Extract the resolution as the last element after splitting the file name
    res <- strsplit(file_name, "_")[[length(strsplit(file_name, "_"))]]  # Get the last element
    res <- res[[length(res)]]  # Get the last element
    res <- gsub(".txt", "", res)  # Remove the file extension
    res <- gsub("R", "", res)  # Remove the "R" prefix

    # Add resolution and tissue type columns
    df$resolution <- res
    df$tissue_type <- tissue_type
    
    # Return the modified data frame
    return(df)
}

#' @name consolidate_data
#' @title Consolidate filtering benchmark results
#' @description This function takes a list of file paths to filtering benchmark result files,
#' reads each file, and consolidates them into a single data frame which is then written to an output file.
#' @param files A character vector containing the paths to the filtering benchmark result files.
#' @return None. The consolidated data frame is written to the output file specified by the global variable OUT.

consolidate_data <- function(files, output_file)  {
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

    fwrite(consolidated_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    message(paste("Consolidated data written to:", output_file))
}