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