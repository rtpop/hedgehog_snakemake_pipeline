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

plot_filtering_bench <- function(df, output_dir, plot_type = "all", plot_title = NULL, plot_file = "filtering_benchmark_plot.pdf") {
  # Load RColorBrewer for color palettes
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please install the 'RColorBrewer' package to use this function.")
  }
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define a color-blind-friendly palette from RColorBrewer
  brewer_palette <- RColorBrewer::brewer.pal(8, "Dark2")
  
  if (plot_type == "all") {
    # Create a facet bar plot for all three metrics
    df_long <- df %>%
      pivot_longer(cols = c("Modularity", "Density", "Number of Edges"),
                   names_to = "Metric",
                   values_to = "Value")
    
    p <- ggplot(df_long, aes(x = factor(resolution), y = Value, fill = Network)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ Metric, scales = "free_x") +  # Allow independent x-axes for each panel
      scale_fill_manual(values = brewer_palette) +  # Apply the RColorBrewer palette
      labs(title = plot_title, x = "Resolution", y = "Value") +
      theme_minimal()
  } else {
    # Create a single bar plot for the specified metric
    p <- ggplot(df, aes_string(x = "factor(resolution)", y = plot_type, fill = "Network")) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = brewer_palette) +  # Apply the RColorBrewer palette
      labs(title = plot_title, x = "Resolution", y = plot_type) +
      theme_minimal()
  }
  
  # Save the plot to a file
  ggsave(plot_file, plot = p, width = 10, height = 6)
}