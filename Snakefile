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
DATA_dir = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
SISANA_DIR = config["sisana_dir"]
SISANA_OUTPUT_DIR = os.path.join(SISANA_DIR, "output") # might want to change later as user can change output dir in the sisana params, but it'll do for now
SRC = config["src_dir"]
ELAND_DIR = os.path.join(OUTPUT_DIR, "eland")
BIHIDEF_DIR = os.path.join(ELAND_DIR, "bihidef")

# Other params
EXP = config["exp_file"]
MOTIF_PRIOR = config["motif_prior_file"]
PPI_PRIOR = config["ppi_prior_file"]
DELIMITER = config["delimiter"]
MAX_COMMUNITIES = config["max_communities"]
MAX_RESOLUTION = config["max_resolution"]

# Input files
SISANA_CONFIG = os.path.join(SISANA_DIR, config["sisana_config"])

# Output files
EXPRESSION_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", EXP + "_filtered.txt")
MOTIF_PRIOR_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", MOTIF_PRIOR + "_filtered.txt")
PPI_PRIOR_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", PPI_PRIOR + "_filtered.txt")
STATS = os.path.join(SISANA_OUTPUT_DIR, "preprocess", EXP + "_filtering_statistics.txt")
PANDA_NET = os.path.join(SISANA_OUTPUT_DIR, "network", "panda_output.txt")
PANDA_EDGELIST = os.path.join(ELAND_DIR, "panda_network_edgelist.txt")
PANDA_NET_FILTERED = os.path.join(ELAND_DIR, "panda_network_filtered.txt")
GENE_COMMUNITIES = os.path.join(BIHIDEF_DIR, "pvg.nodes")


##-------##
## Rules ##
##-------##

## Rule ALL ##
rule all:
    input: 
        #GENE_COMMUNITIES
        PANDA_NET_FILTERED

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

##-------------------------------------##
## Filtering PANDA network for BiHiDef ##
##-------------------------------------##

rule process_and_filter_panda:
    """
    This rule processes and filters the PANDA network.

    Inputs
    ------
    PANDA_NET:
        A TXT file with the PANDA network.
    MOTIF_PRIOR_FILTERED:
        A TXT file with filtered motif prior.
    ------
    Outputs
    -------
    PANDA_EDGELIST:
        A TXT file with the PANDA network processed as an edgelist compatible with networkx.
    PANDA_NET_FILTERED:
        A TXT file with the PANDA edgelist filtered.
    """
    input:
        panda = PANDA_NET, \
        prior = MOTIF_PRIOR_FILTERED
    output:
        edgelist = PANDA_EDGELIST, \
        filtered_net = PANDA_NET_FILTERED
    params:
        process = os.path.join(SRC, "process_panda.py"), \
        filter = os.path.join(SRC, "filter_panda.py"), \
        delimiter = DELIMITER
    container:
        PYTHON_CONTAINER
    message: 
        "; Processing and filtering PANDA network."
    shell:
        """
        echo "Running script with params: script={params.script}, input.panda={input.prior}, output.edgelist={output.edgelist}, delimiter='{params.delimiter}'"
        python {params.process} process_edge_list {input.panda} {output.edgelist} '{params.delimiter}'
        python {params.filter} filter_panda {input.prior} {output.edgelist} {output.filtered_net} '{params.delimiter}'
        """



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
#         script = os.path.join(SRC, "run_bihidef.py"), \
#         max_communities = MAX_COMMUNITIES, \
#         max_resolution = MAX_RESOLUTION
#     message:
#         "; Running BiHiDeF on {input}."
#     shell:
#         """
#         python {params.script} {input} {params.max_communities} {params.max_resolution} {output}
#         """
