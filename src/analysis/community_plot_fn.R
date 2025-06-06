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
  library(data.table)
  library(ggplot2)
  
  df <- fread(data)
  df_unique <- unique(df[, .(tissue, n_selected, n_total)])
  df_unique <- df_unique[order(n_total, decreasing = TRUE)]
  df_unique$tissue <- factor(df_unique$tissue, levels = df_unique$tissue)

  # Reshape to long format for plotting and legend
  df_long <- melt(df_unique, id.vars = "tissue", 
                  measure.vars = c("n_total", "n_selected"),
                  variable.name = "Type", value.name = "Count")
  df_long$Type <- factor(df_long$Type, levels = c("n_total", "n_selected"),
                         labels = c("Total", "Selected"))
  # Set bar width and color for each type
  width_map <- c(Total = 0.7, Selected = 0.5)
  fill_map <- c(Total = "grey80", Selected = "steelblue")
  text_color <- c(Total = "grey30", Selected = "white")

  p <- ggplot(df_long, aes(x = tissue, y = Count, fill = Type, width = width_map[Type])) +
    geom_bar(stat = "identity", position = "identity", color = NA) +
    geom_text(aes(label = Count, color = Type), 
              position = position_identity(), 
              size = 3.5, angle = 90, hjust = 1) +
    scale_fill_manual(values = fill_map, name = "Community Type") +
    scale_color_manual(values = text_color, guide = "none") +
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