#' @name select_top_comm
#' 
#' @param communities String. File path to gmt file with communities.
#' @param scores String. File path for community scores.
#' @param mut_scores String. File path to mutation data. 
#' @param n_comm Integer. Number of "top" communities.

select_top_mut_comm <- function(communities, scores, mut_scores, n_comm = 10) {
    # read data
    comm <- read_gmt(communities)
    comm_scores <- data.table::fread(scores)
    mut <- data.table::fread(mut_scores)

    # select top communities
    top_scores <- apply(comm_scores[,-1], 1, sum)
    names(top_scores) <- comm_scores$V1
    top_scores <- top_scores[order(top_scores, decreasing = TRUE)]
    top_scores <- top_scores[1:n_comm]
    top_comm <- comm[names(top_scores)]

    # get genes mutated in each community
    # top genes
    top_genes <- apply(mut[,-1], 2, sum)
    names(top_genes) <- colnames(mut)[-1]
    mut_genes <- lapply(top_comm, function(x) {
        select_mutated_genes(x, top_genes, n_comm)})
    
    return(mut_genes)
}

#' @name select_mutated_genes
#' @param comm Vector with genes in a community.
#' @param top_genes Vector with sum of mutation scores per gene.
#' @param n_genes Integer. Number of top mutated genes to select.
#' @return List of top mutated genes in each community.
#' 
select_mutated_genes <- function(comm, mut_scores, n_genes = 10) {
    # get in each community
    comm_genes <- mut_scores[comm]
    comm_genes <- comm_genes[order(comm_genes, decreasing = TRUE)]
    comm_genes <- comm_genes[1:n_genes]
    return(names(comm_genes))
}