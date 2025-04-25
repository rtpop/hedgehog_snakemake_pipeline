#' @name plot_n_sig_terms
#' @title Plot the number of significant terms from GO enrichment analysis
#' @description This function takes the results of a GO enrichment analysis and creates a bar plot
#' @param input_file Path to the input file with the results of the GO enrichment analysis
#' @param out_dir Path to output directory
#' @param file_name Name of the output file
#'
#' @return A bar plot showing the number of significant terms for each community 
#' 
plot_n_sig_terms <- function(input_file, out_dir, file_name) {
    
    # Read the input file
    go_res <- data.table::fread(input_file, header = TRUE)

    # Create a heatmap
    p <- ggplot(go_res, aes(x = reorder(community, -n_sig_terms), y = top_go_id, fill = -log10(min_adj_p))) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "blue") +
        theme_minimal() +
        labs(x = "Community", y = "Top GP ID", fill = "-log10(FDR)") +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank())

    # Save the plot
    ggsave(filename = file.path(file_name), plot = p)
}

#' @name plot_top_terms
#' @title Plot the top terms from GO enrichment analysis per community
#' @description This function takes the results of a GO enrichment analysis 
#' and creates a bubble plot with the top term per community.
#' @param input_file Path to the input file with the results of the GO enrichment analysis
#' @param out_dir Path to output directory
#' @param file_name Name of the output file

#' @return A bubble plot showing the top term for each community
#' 
plot_top_terms <- function(input_file, out_dir, file_name) {
    
    go_res <- data.table::fread(input_file, header = TRUE)

    # plot the results
    p <- ggplot(go_res, aes(x = reorder(community, -n_sig_terms), y = term, size = -log10(min_adj_p))) +
        geom_point() +
        theme_minimal() +
        labs(x = "Community", y = "Top Term") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # save the plot
    ggsave(filename = file.path(file_name), plot = p)
}