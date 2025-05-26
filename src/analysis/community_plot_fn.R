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
  # Use only one row per tissue (since n_selected and n_total are repeated)
  df_unique <- unique(df[, .(tissue, n_selected, n_total)])
  # Sort by n_total
  df_unique <- df_unique[order(n_total, decreasing = TRUE)]
  df_unique$tissue <- factor(df_unique$tissue, levels = df_unique$tissue)

  p <- ggplot(df_unique, aes(x = tissue)) +
    geom_bar(aes(y = n_total), stat = "identity", fill = "grey80", width = 0.7) +
    geom_bar(aes(y = n_selected), stat = "identity", fill = "steelblue", width = 0.5) +
    geom_text(aes(y = n_total, label = n_total), hjust = 1, color = "grey30", size = 3.5, angle = 90) +
    geom_text(aes(y = n_selected, label = n_selected), hjust = 1, color = "white", size = 3.5, angle = 90) +
    labs(title = "Communities per Tissue",
         x = "Tissue",
         y = "Number of Communities") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.direction = "horizontal")

  if (is.null(file_name)) {
    return(p)
  } else {
    ggsave(file_name, plot = p, width = 8, height = 5)
  }
}