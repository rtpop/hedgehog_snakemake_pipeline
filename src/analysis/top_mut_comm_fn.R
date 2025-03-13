#' @name select_top_comm
#' 
#' @param communities String. File path to gmt file with communities.
#' @param scores String. File path for community scores.
#' @param mut String. File path to mutation data. 
#' @param n_comm Integer. Number of "top" communities.

select_top_mut_comm <- function(communities, scores, mut, n_comm = 10) {
    # read data
    comm <- read_gmt(communities)
    comm_scores <- data.table::fread(scores)
    mutation <- data.table::fread(mut)
    samples <- mutation$V1
    mutation <- mutation[, -1]
    mutation <- as.matrix(mutation)
    rownames(mutation) <- samples

    # select top communities
    top_scores <- apply(comm_scores[,-1], 1, mean)
    names(top_scores) <- comm_scores$V1
    top_scores <- top_scores[order(top_scores, decreasing = TRUE)]
    top_scores <- top_scores[1:n_comm]
    top_comm <- comm[names(top_scores)]

    # get patients that have mutations in the top communities
    patients_with_mut <- lapply(names(top_comm), function(comm_name) {
        comm_index <- which(comm_scores$V1 == comm_name)
        comm_scores_subset <- comm_scores[comm_index, -1, with = FALSE]
        patients <- colnames(comm_scores_subset)[apply(comm_scores_subset, 2, function(x) any(x > 0))]
        return(patients)
    })
    names(patients_with_mut) <- names(top_comm)

    # get mutated genes for each patient in each community
    mut_genes <- lapply(names(top_comm), function(comm_name) {
        comm_genes <- top_comm[[comm_name]]
        comm_patients <- patients_with_mut[[comm_name]]
        select_mutated_genes(comm_genes, mutation, comm_patients)
    })
    names(mut_genes) <- names(top_comm)
    
    return(mut_genes)
}

#' @name select_mutated_genes
#' @param comm_genes Vector with genes in a community.
#' @param mutation Data frame with mutation data.
#' @param patients Vector with patient IDs.
#' @return List of mutated genes for each patient in the community.
#' 
select_mutated_genes <- function(comm_genes, mutation, patients) {
    # Subset mutation data for the given community and patients
    idx <- intersect(colnames(mutation), comm_genes)
    comm_mutation <- mutation[patients, idx, drop = FALSE]
    
    # Find mutated genes for each patient
    mutated_genes <- lapply(patients, function(patient) {
        genes <- comm_genes[comm_mutation[patient, ] > 0]
        return(genes)
    })
    names(mutated_genes) <- patients
    
    return(mutated_genes)
}

#' @name mut_gene_summary
#' @param mut_genes List of mutated genes for each patient in the community.
#' @return Data frame with summary of mutated genes.
#'
mut_gene_summary <- function(mut_genes) {
    # Get unique genes
    all_genes <- unique(unlist(mut_genes))
    
    # Initialize gene summary data frame
    gene_summary <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
    
    # Count number of mutations for each gene
    gene_counts <- lapply(mut_genes, function(comm) {
        count <- table(unlist(comm))
        return(count)
    })
    
    # Convert gene_counts list to data frames
    gene_counts_dfs <- lapply(seq_along(gene_counts), function(i) {
        df <- data.frame(gene = names(gene_counts[[i]]), count = as.numeric(gene_counts[[i]]), stringsAsFactors = FALSE)
        colnames(df)[2] <- names(gene_counts[i])
        return(df)
    })
    
    # Perform a full join on all data frames in gene_counts_dfs
    gene_counts_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "gene"), gene_counts_dfs)
    
    # Merge gene_summary with gene_counts_df
    gene_summary <- dplyr::full_join(gene_summary, gene_counts_df, by = "gene")
    
    return(gene_summary)
}