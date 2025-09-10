
#' @name get_gtex_data
#' 
#' @param data RData file containing GTEx data.
#' @param extract_edges Logical. If edges should be extracted for each tissue. If FALSE,
#' expression will be extracted instead

get_gtex_data <- function(data, extract_edges = TRUE, out_dir = "data/") {
    # Load the data
    load(data)

    tissues <- as.character(unique(samples$Tissue))
    # Get the gene names
    gene_names <- genes$Symbol[match(edges[,2], genes$Name)]
    edges[,2] <- gene_names
    data.table::fwrite(edges, file = file.path(out_dir, "motif_prior.txt"), sep = ",", row.names = FALSE, col.names = FALSE)

    for (i in tissues) {
        tissue_samples <- samples[samples$Tissue == i, 1]
        

        if (extract_edges) {
            # Extract edges for the specific tissue
            res <- net[, i]
            res <- cbind(edges[,1:2], res)
            file_name <- paste0("panda_network_edgelist.txt")
        } else {
            # Extract expression data for the specific tissue
            res <- data.frame(exp[, tissue_samples])
            gene_names <- as.character(genes[which(genes$Name %in% res$Gene), 1])
            res <- cbind(gene_names, res)
            file_name <- paste0("exp_", i, ".txt")
        }

        # Create directory if it doesn't exist
        dir_path <- file.path(out_dir, i)
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }

        # Save the data with row names
        data.table::fwrite(res, file = file.path(dir_path, file_name), sep = ",", row.names = FALSE, col.names = FALSE)
    }
}