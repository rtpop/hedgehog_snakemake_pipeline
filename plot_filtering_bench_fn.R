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
  
  # Add resolution and tissue type columns
    res <- strsplit(file_name, "_")[[1]][[3]]
    res <- gsub(".txt", "", res)
    res <- gsub("R", "", res)

  df$resolution <- res
  df$tissue_type <- tissue_type
  
  # Return the modified data frame
  return(df)
}


#' @name plot_filtering_bench
#' 
#' @title Plot filtering benchmark results
#' 
#' @description This function takes a data frame containing filtering benchmark results and generates a plot comparing 
#' the performance of different filtering methods.
#' 
#' @param df A data frame containing the filtering benchmark results. The data frame should have the following columns:
#' \itemize{
#'  \item \code{Network}: The method of filtering used (ELAND, Prior or Top).
#'  \item \code{Modularity}: The modularity score of the filtered network.
#'  \item \code{Density}: The density of the filtered network.
#'  \item \code{Number of Edges}: The number of edges in the filtered network.
#'  \item \code{resolution}: The resolution used in the filtering benchmark.
#'  \item \code{tissue_type}: The type of tissue used in the filtering benchmark.
#' }
#' @param output_dir The directory where the plot will be saved.
#' @param plot_type The type of plot to generate. This should be one of "modularity", "density", or "edges".
#' @param plot_title The title of the plot.
#' 

plot_filtering_bench <- function(df, output_dir, plot_type, plot_title=NULL) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set the file name for the plot
  file_name <- paste0(output_dir, "/", plot_type, "_", gsub(" ", "_", plot_title), ".pdf")
  
  # Create the plot
  p <- ggplot(df, aes_string(x = "resolution", y = plot_type, color = "Network")) +
    geom_point() +
    geom_line() +
    labs(title = plot_title, x = "Resolution", y = plot_type) +
    theme_minimal()
  
  # Save the plot to a file
  ggsave(file_name, plot = p)
}