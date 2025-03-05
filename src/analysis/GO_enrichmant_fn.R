#' @name Read GMT file
#' 
#' @param gmt_file Path to the GMT file
#' 
#' @return A list with the gene sets
#' 
#' @examples
#' gene_sets <- read_gmt("path/to/gmt_file.gmt")

read_gmt <- function(gmt_file) {
    # Read the contents of the GMT file into a character vector
    gmt <- readLines(gmt_file)
    
    # Initialize an empty list to store the gene sets
    gene_sets <- list()
    
    # Loop through each line in the GMT file
    for (line in gmt) {
        # Split the line by tab characters
        line_split <- strsplit(line, "\t")[[1]]
        
        # Remove the second element (description) from the split line
        line_split <- line_split[-2]
        
        # Use the first element as the gene set name and the remaining elements as the genes
        gene_sets[[line_split[1]]] <- line_split[-1]
    }
    
    # Return the list of gene sets
    return(gene_sets)
}

#' @description RUN GO term enrichment analysis on BiHiDef communities
#' 
#' @param gmt_file GMT file with the communities.
#' @param out_dir Path to output directory
#' @param auto_bg Logical, if TRUE, the background is automatically set to the genes in the gene sets.
#' @param bg_file Optional. Path to a file with the background genes. Only needed if auto_bg is FALSE.
#' @param save_all Logical, if TRUE the GO enrichment for each community will be saved in addition to 
#' a summary for all communities. Otherwise, only the summary will be saved.
#' @param sig_thresh Threshold for significance of enrichment. This is applied to the adjusted p value.
#' @param statistic String. What test statistic should be used. See \code{\link{topGO}} for available options.
#' @param algorithm String. What algorithm should be used for enrichment. See \code{\link{topGO}} for available options.
#'
#' @return A data frame with the results of the GO term enrichment analysis.
#'
#' @examples
#' gene_list <- c("gene1", "gene2", "gene3")
#' gene_sets <- list("gene_set1" = c("gene1", "gene2"), "gene_set2" = c("gene2", "gene3"))
#' go_results <- run_go_enrichment(gene_sets, gene_list)
#' print(go_results)


run_go_enrichment <- function(gmt_file, out_dir, auto_bg = TRUE, bg_file = NULL, save_all = TRUE, sig_thresh = 0.05,
                            statistic = "fisher", LGORITHM = "classic") {
    # Read the gmt file
    gene_sets <- read_gmt(gmt_file)

    # get all the genes in the gene sets
    gene_list <- unique(unlist(gene_sets))

    # set background genes
    if (auto_bg) {
        bg_genes <- gene_list
    } else {
        anno <- data.table::fread(bg_file, header = TRUE)
        bg_genes <- anno$symbol
    }

    # mapping GO terms to genes
    go_gene_mapping <- topGO::annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
	gene_go_mapping <- topGO::inverseList(go_gene_mapping)

    # run enrichment for each community
    res_list <- lapply(names(gene_sets), function(comm_name) {
        gene_set <- gene_sets[[comm_name]]
        enrich_community(gene_set = gene_set, background = bg_genes, mapping = gene_go_mapping, 
                        save_all = save_all, comm_name = comm_name, sig_thresh = sig_thresh,
                        algorithm = algorithm, statistic = statistic)
    })

    #Combine results into a single data frame
    res_df <- do.call(rbind, lapply(res_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
    colnames(res_df) <- c("community", "n_sig_terms", "min_adj_p", "top_go_id", "top_go_term")
    
    # write file
    data.table::fwrite(res_df, sep="\t", file = file.path(out_dir, paste0(statistic, "_GO_Summary.txt")))
}

#' @name enrich_community
#' 
#' @inheritParams run_go_enrichment
#' @param gene_set
#' @param background
#' @param mapping
#' @param comm_name Name of community being tested
#' 

enrich_community <- function(gene_set, background, mapping, save_all, comm_name, out_dir, sig_thresh, statistic, algorithm) {
    # get the genes from the gene set that are in the bg
    # I feel like this is somewhat redundant because I don't think it's possible to
    # have genes in the communities that are not in the bg
    # but this is how it was done in M's og script, so for now...
    interesting <- background[which(background %in% gene_set)]

    # create the allGenes and geneSel vectors
    # in this case they are the same, so I will only make one of them
    selection = factor(as.integer(background %in% interesting))
    names(selection) <- background

    # create topGO object
    go_object <- topGO::new("topGOdata", ontology = "BP", allGenes = selection, geneSel = selection,
                            annot = annFUN.gene2GO, gene2GO = mapping)

    # run fisher test
    go_fisher <- topGO::runTest(go_object, algorithm = algorithm, statistic = statistic)

    # generate results table
    enrich_res <- topGO::GenTable(go_object, classic = go_fisher, orderBy = "classic", ranksOf = "classic")
    p_adj <- p.adjust(enrich_res$classic, "BH")
    enrich_res$p_adj <- p_adj
    enrich_res <- enrich_res[order(enrich_res)]

    if (save_all) {
        data.table::fwrite(enrich_res, file = file.path(out_dir, paste0(comm_name, ".txt")),
                            row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
    }

    # getting the summarised results
    n_sig_term <- length(which(enrich_res$p_adj < sig_thresh))
    min_adj_p <- min(enrich_res$p_adj)
    top_go_id <- enrich_res[1, 1]
    top_go_term <- enrich_res[1, 2]
    res <- c(comm_name, n_sig_term, min_adj_p, top_go_id, top_go_term)

    return(res)
}