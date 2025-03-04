##--------------------------------------------------------------------------##
## HOW TO RUN                                                               ##
## Run from directory containing Snakefile                                  ##
##--------------------------------------------------------------------------##
## For dry run                                                              ##
## snakemake --cores 1 -np                                                  ##
##--------------------------------------------------------------------------##
## For local run                                                            ##
## snakemake --cores 1                                                      ##
##--------------------------------------------------------------------------##
## For running with singularity container                                   ##
## snakemake --cores 1 --use-singularity --singularity-args '\-e' --cores 1 ##
##--------------------------------------------------------------------------##

##-----------##
## Libraries ##
##-----------##

import os 
import sys
import glob
from pathlib import Path
import time

##------------##
## Parameters ##
##------------##

# Config file
global CONFIG_PATH
CONFIG_PATH = "config.yaml"
configfile: CONFIG_PATH

# Containers
PYTHON_CONTAINER = config["python_container"]
ANALYSIS_CONTAINER = config["analysis_container"]

# Directories
DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
SISANA_DIR = config["sisana_dir"]
SISANA_OUTPUT_DIR = os.path.join(SISANA_DIR, "output") # might want to change later as user can change output dir in the sisana params, but it'll do for now
SRC = config["src_dir"]
ELAND_DIR = os.path.join(OUTPUT_DIR, "eland")
BIHIDEF_DIR = os.path.join(ELAND_DIR, "bihidef")
SAMBAR_DIR = os.path.join(ELAND_DIR, "sambar")


# BiHiDeF params
TAR_TAG = config["target_tag"]
REG_TAG = config["regulator_tag"]
MAX_COMMUNITIES = config["max_communities"]
MAX_RESOLUTION = config["max_resolution"]
MAX_GENES = config["max_genes"]
MIN_GENES = config["min_genes"]
BIHIDEF_RUN_DIR = os.path.join(BIHIDEF_DIR, "C" + str(MAX_COMMUNITIES) + "_R" + str(MAX_RESOLUTION)) # for running with different params

# Other params
EXP = config["exp_file"]
MOTIF_PRIOR = config["motif_prior_file"]
PPI_PRIOR = config["ppi_prior_file"]
DELIMITER = config["delimiter"]

# Input files
SISANA_CONFIG = os.path.join(SISANA_DIR, config["sisana_config"])
MUT_DATA = os.path.join(DATA_DIR, config["mut_data"])
ESIZE = os.path.join(DATA_DIR, config["esize_file"])
CAN_GENES = os.path.join(DATA_DIR, config["can_genes"])

# samba params
SAMBAR_OUTPUT_DIR = os.path.join(SAMBAR_DIR, "C" + str(MAX_COMMUNITIES) + "_R" + str(MAX_RESOLUTION))

# sisana outputs
EXPRESSION_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", EXP + "_filtered.txt")
MOTIF_PRIOR_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", MOTIF_PRIOR + "_filtered.txt")
PPI_PRIOR_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", PPI_PRIOR + "_filtered.txt")
STATS = os.path.join(SISANA_OUTPUT_DIR, "preprocess", EXP + "_filtering_statistics.txt")
PANDA_NET = os.path.join(SISANA_OUTPUT_DIR, "network", "panda_network.txt")

# panda processing outputs
PANDA_EDGELIST = os.path.join(ELAND_DIR, "panda_network_edgelist.txt")
PANDA_NET_FILTERED = os.path.join(ELAND_DIR, "panda_network_filtered.txt")

# BiHiDeF outputs
GENE_COMMUNITIES = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG  + ".nodes")
SELECTED_COMMUNITIES = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + "_selected_communities.gmt")
COMMUNITY_STATS = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + "_community_stats.txt")

# sambar outputs
PATHWAY_SCORES = os.path.join(SAMBAR_DIR, "pt_out.csv")
MUTATION_SCORES = os.path.join(SAMBAR_DIR, "mt_out.csv")

# Downstream analysis outputs
GO_ENRICHMENT = os.path.join(ELAND_DIR, "go_enrichment.txt")

##-------##
## Rules ##
##-------##

## Rule ALL ##
rule all:
    input: 
        SELECTED_COMMUNITIES, \
        #PANDA_NET_FILTERED, \
        PATHWAY_SCORES


##-----------------------##
## Making PANDA networks ##
##-----------------------##

# rule run_sisana:
#     """
#     This rule runs PANDA using the SiSaNA pipeline.

#     SiSaNA is available at
#     https://github.com/kuijjerlab/sisana

#     Inputs
#     ------
#     SISANA_CONFIG:
#         Config yml file for SiSaNA.
#     ------
#     Outputs
#     -------
#     EXPRESSION_FILTERED:
#         A TXT file with filtered expression.
#     MOTIF_PRIOR_FILTERED:
#         A TXT file with filtered motif prior.
#     PPI_PRIOR_FILTERED:
#         A TXT file with filtered PPI prior.
#     STATS:
#         A file containing information about genes filtered.
#     PANDA_NET:
#         A TXT file with the PANDA network.
#     """
#     input:
#         SISANA_CONFIG
#     output:
#         EXPRESSION_FILTERED, \
#         MOTIF_PRIOR_FILTERED, \
#         PPI_PRIOR_FILTERED, \
#         STATS, \
#         PANDA_NET
#     container:
#         PYTHON_CONTAINER
#     message: 
#         "; Running sisana preprocess on {input}."
#     shell:
#         """
#         sisana preprocess {input}
#         sisana generate {input}
#         """

# ##-------------------------------------##
# ## Filtering PANDA network for BiHiDef ##
# ##-------------------------------------##

# rule process_and_filter_panda:
#     """
#     This rule processes and filters the PANDA network.

#     Inputs
#     ------
#     PANDA_NET:
#         A TXT file with the PANDA network.
#     MOTIF_PRIOR_FILTERED:
#         A TXT file with filtered motif prior.
#     ------
#     Outputs
#     -------
#     PANDA_EDGELIST:
#         A TXT file with the PANDA network processed as an edgelist compatible with networkx.
#     PANDA_NET_FILTERED:
#         A TXT file with the PANDA edgelist filtered.
#     """
#     input:
#         panda = PANDA_NET, \
#         prior = MOTIF_PRIOR_FILTERED
#     output:
#         edgelist = PANDA_EDGELIST, \
#         filtered_net = PANDA_NET_FILTERED, \
#         updated_prior_sep = MOTIF_PRIOR_FILTERED.replace(".txt", "_updated_sep.txt"), \
#         updated_panda_sep = PANDA_NET.replace(".txt", "_updated_sep.txt") 
#     params:
#         process = os.path.join(SRC, "process_panda.py"), \
#         filter = os.path.join(SRC, "filter_panda.py"), \
#         delimiter = DELIMITER
#     container:
#         PYTHON_CONTAINER
#     message: 
#         "; Processing and filtering PANDA network." \
#         "Running {params.process} on {input.panda} and {input.prior} to create {output.edgelist} with --delimiter {params.delimiter}." \
#         "Running {params.filter} on {output.updated_panda_sep} and {input.prior} to create {output.filtered_net} with --delimiter {params.delimiter}."
#     shell:
#         """
#         python {params.process} {input.panda} {input.prior} {output.edgelist} --delimiter '{params.delimiter}'
#         head {output.updated_prior_sep}
#         python {params.filter} {output.updated_prior_sep} {output.edgelist} {output.filtered_net} --delimiter '{params.delimiter}'
#         """

# ## --------------- ##
# ## Running BiHiDeF ##
# ## --------------- ##

# rule run_bihidef:
#     """
#     This rule runs the BiHiDeF algorithm.

#     BiHiDeF is available at
#     """
#     input:
#         PANDA_NET_FILTERED
#     output:
#         gene_communities = GENE_COMMUNITIES
#     params:
#         run_script = os.path.join(SRC, "run_bihidef.py"), \
#         measure_script = os.path.join(SRC, "measure_resources.py"), \
#         max_communities = MAX_COMMUNITIES, \
#         max_resolution = MAX_RESOLUTION, \
#         output_prefix_reg = REG_TAG, \
#         output_prefix_target = TAR_TAG, \
#         out_dir = BIHIDEF_RUN_DIR, \
#         log_file = os.path.join(BIHIDEF_RUN_DIR, "run_log.log")
#     container:
#         PYTHON_CONTAINER
#     message:
#         "; Running BiHiDeF on {input} with params:" \
#             "--comm_mult {params.max_communities}" \
#             "--max_res {params.max_resolution}" \
#             "--output_dir {params.out_dir}" \
#             "--output_prefix_reg {params.output_prefix_reg}" \
#             "--output_prefix_tar {params.output_prefix_target}"
#     shell:
#         """
#         mkdir -p {params.out_dir}
#         python {params.measure_script} {params.log_file} "python {params.run_script} {input} --comm_mult {params.max_communities} --max_res {params.max_resolution} \
#         --output_dir {params.out_dir} --output_prefix_reg {params.output_prefix_reg} --output_prefix_tar {params.output_prefix_target}"
#        """

## --------------------- ##
## Selecting communities ##
## --------------------- ##

rule select_communities:
    """
    This rule selects the communities from the BiHiDeF output and formats them as a GMT file.

    Inputs
    ------
    GENE_COMMUNITIES:
        A TXT file with the communities from BiHiDeF.
    ------
    Outputs
    -------
    SELECTED_COMMUNITIES:
        A TXT file with the selected communities.
    COMMUNITY_STATS:
        A TXT file with statistics about the communities.
    """
    input:
        GENE_COMMUNITIES
    output:
        selected_communities = SELECTED_COMMUNITIES, \
        stats = COMMUNITY_STATS
    params:
        script = os.path.join(SRC, "select_communities.py"), \
        max_genes = MAX_GENES, \
        min_genes = MIN_GENES
    container:
        PYTHON_CONTAINER
    message:
        "; Selecting communities from {input} with params:" \
            "--max_size {params.max_genes}" \
            "--min_size {params.min_genes}" \
            "--log {output.stats}"
            "output {output.selected_communities}"
    shell:
        """
        python {params.script} {input} {output.selected_communities} --log {output.stats} --max_size {params.max_genes} --min_size {params.min_genes}
        """

## -------------- ##
## Running sambar ##
## -------------- ##

rule run_sambar:
    """
    This rule runs the sambar algorithm.

    Sambar is available at:
    """
    input:
        gmt = SELECTED_COMMUNITIES, \
        mut = MUT_DATA, \
        esize = ESIZE, \
        can_genes = CAN_GENES
    output:
        pathway_scores = PATHWAY_SCORES, \
        mutation_scores = MUTATION_SCORES

    params:
        script = os.path.join(SRC, "run_sambar.py"), \
        out_dir = SAMBAR_OUTPUT_DIR
    container:
        PYTHON_CONTAINER
    message:
        "; Running sambar on {input} with params:" \
            "--output-dir {params.out_dir} " \
            "--gmt-file {input.gmt} " \
            "--mut-file {input.mut} " \
            "--esize-file {input.esize} " \
            "--can-genes {input.can_genes} "
    shell:
        """
        mkdir -p {params.out_dir}
        python {params.script} --gmt-file {input.gmt} --mut-file {input.mut} --esize-file {input.esize} --output-dir {params.out_dir} --can-genes {input.can_genes}
        """

## ---------------------------- ##
## GO enrichment of communities ##
## ---------------------------- ##

rule go_enrichment:
    """
    This rule runs GO enrichment on the selected communities.

    Inputs
    ------
    SELECTED_COMMUNITIES:
        A GMT file with the selected communities.
    ------
    Outputs
    -------
    GO_ENRICHMENT:
        A TXT file with the GO enrichment results.
    """
    input:
        SELECTED_COMMUNITIES
    output:
        go_enrichment = os.path.join(ELAND_DIR, "go_enrichment.txt")
    params:
        script = os.path.join(SRC, "run_go_enrichment.py")
    container:
        ANALYSIS_CONTAINER
    message:
        "; Running GO enrichment on {input} with params:" \
            "--output {output.go_enrichment}"
    shell:
        """
        Rscript {params.script} {input} {output.go_enrichment}
        """