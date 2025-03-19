## ----------------------------------------------------------------------------------- ##
## This script is not run as part of the snakemake and it is specific to the GTEx data ##
## Found on Zenodo https://zenodo.org/records/838734                                   ##
## ----------------------------------------------------------------------------------- ##

#' @name get_gtex_exp
#' 
#' @param data RData object with gtex tissues. Downloaded from Zenodo.
#' 
get_gtex_exp <- function(data){
    load(data)

    tissues <- as.character(unique(samples$Tissue))

    for (i in tissues) {
        tissue_samples <- samples[samples$Tissue == i,1]
        tissue_exp <- exp[,tissue_samples]
        gene_names <- genes[rownames(tissue_exp),1]
        tissue_exp <- tissue_exp[order(rownames(tissue_exp)),]

        data.table::fwrite(tissue_exp, file = file.path("data", i, paste0("gtex_", i, "_exp.tsv"), sep = "\t"))
    }
}

data <- "/storage/kuijjerarea/romana/eland/ELAND/data/GTEx_PANDA_net.RData"
get_gtex_exp(data)