#' @name process_community
#' @title Process community data
#' @description Process community data
#' @param gmt_file Path to GMT file
#' @param stats_file Path to stats file
#' @param tissue Tissue name
#' 
#' @return Data frame with processed community data

# process_community <- function(gmt_file, stats_file, tissue) {
#   # Read GMT file
#   gmt_data <- read_gmt(gmt_file)

#   # Read stats file
#   stats <- readLines(stats_file)

#     # Skip the first line (time and date)
#   stats <- stats[-1]

#   # Extract numbers from each line
#   numbers <- as.numeric(gsub(".*: ", "", stats))

#   # Create a data frame with the community stats
#     df <- data.frame(
#         tissue = tissue,
#         n_comm = numbers[3],
#         n_comm_total = numbers[4],
#         max_size = numbers[1],
#         min_size = numbers[2],
#         tissue
#     )

#     comm_size <- data.frame(
#         comm = names(gmt_data),
#         size = sapply(gmt_data, length)
#     )

#     return(list(df = df, comm_size = comm_size))
# }

#' @name plot_n_communities
#' @title Plot number of communities
#' @description Plot number of communities
#' @param data Data frame containing community stats across tissues
#' @param file_name Path to output file
#' @return Plot object or NULL

plot_n_communities <- function(data file_name = NULL) {
  
  # Plot number of selected communities for all tissues
  p <- ggplot(df, aes(x = tissue, y = n_selectd)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Number of Selected Communities per Tissue",
         x = "Tissue",
         y = "Number of Selected Communities") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (is.null(file_name)) {
    return(p)
  } else {
    ggsave(file_name, plot = p, width = 8, height = 5)
  }
}

#' @name plot_community_size
#' @title Plot community size
#' @description Plot community size
#' @param data Data frame containing the community size
#' @param tissue Tissue name
#' @param file_name Path to output file
#' @return Plot object or NULL

plot_community_size <- function(data, tissue, file_name = NULL) {
  # Create a bar plot of the community size
  p <- ggplot(data, aes(x = tissue, y = size)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Community Size in", tissue),
         x = "Tissue",
         y = "Community Size") +
    theme_minimal()

  if (is.null(file_name)) {
    # If no file name is provided, display the plot
    return(p)
  } else {
    # Save the plot to a file
    ggsave(file_name, plot = p)
  }
}