#' @name cluster_scores
#' @param scores Character. Path to file with community/pathway scores as output by PySambar.
#' @param out_dir Character. Parth to output directory.
#' @param distance_matrix Character. Path to file with distance matrix as output by PySambar.
#' @param binaries Logical. If TRUE, the pathway scores will be binarised for plotting.
#' The clustering will always be done on the un-modified scores, this is just for visualisation.
#' @param log_transform Logical. If TRUE, the pathway scores will be -log10 transformed for plotting.
#' The clustering will always be done on the un-modified scores, this is just for visualisation.
#' 

cluster_scores <- function(scores, out_dir, distance_matrix, binarise = TRUE, log_transform = FALSE) {
    # sanity check
    if (binarise == TRUE && log_transform == TRUE) {
        stop("binarise and log_transform cannot both be TRUE. Please set one or both to FALSE.")
    }
    
    # read scores
    comm_scores <- data.table::fread(scores)

    # read distance matrix
    samples <- colnames(comm_scores)[-1]
    dist <- prepare_distance_matrix(distance_matrix, samples)

    # convert to matrix
    communities <- comm_scores[, 1]
    comm_mat <- as.matrix(comm_scores[, -1])
    rownames(comm_mat) <- NULL
    colnames(comm_mat) <- NULL

    # cluster
    comm_clust <- stats::hclust(comm_dist)

    # transform scores
    if (binarise) {
        comm_mat <- ifelse(comm_mat > 0, 1, 0)
        tag <- "bin"
        # set palette
        bw_palette <- colorRampPalette(c("white", "black"))(2)
    } else if (log_transform) {
        comm_mat <- -log10(comm_mat + 1)
        bw_palette <- colorRampPalette(c("white", "black"))(100)
        tag <- "log"
    } else {
        tag <- "raw"
    }

    # plot heatmap
    comm_heat <- pheatmap::pheatmap(comm_mat, cluster_rows = TRUE, cluster_cols = comm_clust, labels_row = NULL, labels_col = NULL, color = bw_palette)
    grid::grid.text("samples", x = unit(0.5, "npc"), y = unit(-0.03, "npc"), gp = gpar(fontsize = 12))
    grid::grid.text("communities", x = unit(-0.03, "npc"), y = unit(0.5, "npc"), rot = 90, gp = gpar(fontsize = 12))
    ggplot2::ggsave(comm_heat, file = file.path(out_dir, "community_scores_clusters_", tag, "_.pdf"))
}

#' @name prepare_distance_matrix
#' @inheritParams cluster_scores
#' @param samples Character vector with sample names. 

prepare_distance_matrix <- function(distance_matrix, samples) {
    # read file
    dist <- data.table::fread(distance_matrix)

    # set names
    colnames(dist) <- communities
    rownames(dist) <- communities

    # convert to dist object
    dist <- as.matrix(dist)
    dist <- stats::as.dist(dist)

    return(dist)
}