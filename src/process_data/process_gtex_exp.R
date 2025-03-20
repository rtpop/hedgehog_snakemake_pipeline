## ----------------------------------------------------------------------------------- ##
## This script is not run as part of the snakemake and it is specific to the GTEx data ##
## Found on Zenodo https://zenodo.org/records/838734                                   ##
## ----------------------------------------------------------------------------------- ##

rm(list=ls())

#' @name get_gtex_exp
#' 
#' @param exp Data frame with gene expression data.
#' @param samples Data frame with sample information.
#' @param genes Data frame with gene information.
#' 
get_gtex_exp <- function(exp, samples, genes) {

    tissues <- as.character(unique(samples$Tissue))

    for (i in tissues) {
        tissue_samples <- samples[samples$Tissue == i, 1]
        tissue_exp <- data.frame(exp[, tissue_samples])
        gene_names <- as.character(genes[which(genes$Name == rownames(tissue_exp)), 1])
        tissue_exp <- cbind(gene_names, tissue_exp)

        # Create directory if it doesn't exist
        dir_path <- file.path("data", i)
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }

        # Save the data with row names
        data.table::fwrite(tissue_exp, file = file.path(dir_path, paste0("gtex_", i, "_exp.tsv")), sep = "\t")
    }
}

data <- "/storage/kuijjerarea/romana/eland/ELAND/data/GTEx_PANDA_net.RData"
load(data)
get_gtex_exp(exp, samples, genes)