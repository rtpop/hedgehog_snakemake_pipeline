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

plot_n_communities <- function(data, file_name = NULL) {
  df <- fread(data)

  # Ensure n_selected is numeric
  df$n_selected <- as.numeric(df$n_selected)
  # Sort by n_selected
  df <- df[order(df$n_selected), ]
  # Set tissue as a factor in the sorted order
  df$tissue <- factor(df$tissue, levels = unique(df$tissue))

  # Plot number of selected communities for all tissues
  p <- ggplot(df, aes(x = tissue, y = n_selected)) +
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