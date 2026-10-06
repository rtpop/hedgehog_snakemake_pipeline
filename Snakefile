## ------------------------------------------------------------------------------------------- ##
## HOW TO RUN                                                                                  ##
## Run from directory containing Snakefile                                                     ##
## ------------------------------------------------------------------------------------------- ##
## For dry run                                                                                 ##
## snakemake --cores 1 -np                                                                     ##
## ------------------------------------------------------------------------------------------- ##
## For local run                                                                               ##
## snakemake --cores 1                                                                         ##
## ------------------------------------------------------------------------------------------- ##
## For running with singularity container                                                      ##
## snakemake --cores 1 --use-singularity --singularity-args '\-e' --cores 1                    ##
## ------------------------------------------------------------------------------------------- ## ------------------------------ ##
## For running in the background                                                                                                 ##
## nohup snakemake --cores 1 --use-singularity --singularity-args '\-e' 2>&1 > logs/snakemake_$(date +'%Y-%m-%d_%H-%M-%S').log & ##
## ----------------------------------------------------------------------------------------------------------------------------- ##

##-----------##
## Libraries ##
##-----------##

import os 
import sys
import glob
from pathlib import Path
import time

## ----------------- ##
## Global parameters ##
## ----------------- ##

configfile: "config.yaml"

NCORES = config["ncores"]

# Containers
PYTHON_CONTAINER = config["python_container"]
R_CONTAINER = config["r_container"]

# Directories
DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
SRC = config["src_dir"]
HEDGEHOG_DIR = os.path.join(OUTPUT_DIR, "{tissue_type}", "hedgehog_bug_fix_q_score")
BIHIDEF_RUN_DIR = os.path.join(HEDGEHOG_DIR, "bihidef_run")
GO_DIR = os.path.join(OUTPUT_DIR, "{tissue_type}", "go_enrichment")

# Other params
DELIMITER = config["delimiter"]
TISSUE = config["tissue"] # wildcard
TAR_TAG = config["target_tag"]

## ------------------------------ ##
## Download and process GTEx data ##
## ------------------------------ ##
GTEX_DATA_FILE = os.path.join(DATA_DIR, "download", config["gtex_data_file"])
PROCESSING_LOG = os.path.join(DATA_DIR, config["processing_log"])
MOTIF_PRIOR = os.path.join(DATA_DIR, config["motif_prior"])

##-------------------------------------##
## Filtering PANDA network for BiHiDef ##
##-------------------------------------##
PANDA_NET = os.path.join(DATA_DIR, "{tissue_type}", config["panda_net_file"])
FILTERING_METHOD = config["filtering_method"]


## ----------- ##
## Run BiHiDef ##
## ----------- ##
GENE_COMMUNITIES = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + ".nodes")
MAX_COMMUNITIES = config["max_communities"]
MAX_RESOLUTION = config["max_resolution"]
FILTERING_BIHIDEF = config["filtering_bihidef"]
REG_TAG = config["regulator_tag"]
TAR_TAG = config["target_tag"]


# set filtering params
if FILTERING_METHOD == "both":
    PANDA_NET_FILTERED = [
        os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt"),
        os.path.join(HEDGEHOG_DIR, "panda_network_filtered_hedgehog.txt")
    ]
elif FILTERING_METHOD == "prior":
    PANDA_NET_FILTERED = [os.path.join(HEDGEHOG_DIR, "panda_network_filtered_prior.txt")]
elif FILTERING_METHOD == "hedgehog":
    PANDA_NET_FILTERED = [os.path.join(HEDGEHOG_DIR, "panda_network_filtered_hedgehog.txt")]
else:
    raise ValueError("Unknown filtering method: {}".format(FILTERING_METHOD))

##-------##
## RULES ##
##-------##

rule all:
    input:
        expand(PANDA_NET_FILTERED, tissue_type=TISSUE),
        expand(GENE_COMMUNITIES, tissue_type=TISSUE)
## ---------------------------- ##
## Download & process GTEX data ##
## ---------------------------- ##

rule downaload_gtex:
    output:
        GTEX_DATA_FILE
    params:
        out_dir = os.path.join(DATA_DIR, "download"), \
        log_file = os.path.join(DATA_DIR, "download", "download_gtex.log")
    container:
        PYTHON_CONTAINER
    message:
        "; Downloading GTEx data."
    shell:
        """
        mkdir -p {params.out_dir}
        curl --output {output} https://zenodo.org/records/838734/files/GTEx_PANDA_tissues.RData?download=1 > {params.log_file} 2>&1
        """

rule process_gtex:
    input:
        gtex_data = GTEX_DATA_FILE
    output:
        output_log = PROCESSING_LOG,
        prior = MOTIF_PRIOR
    params:
        out_dir = DATA_DIR,
        extract_edges = True
    container:
        R_CONTAINER
    message:
        "; Processing GTEx data."
    script:
        os.path.join(SRC, "process_gtex.R")

##-------------------------------------##
## Filtering PANDA network for BiHiDef ##
##-------------------------------------##
rule process_and_filter_panda:
    input:
        panda = PANDA_NET,
        prior = MOTIF_PRIOR
    output:
        filtered_panda = PANDA_NET_FILTERED
    params:
        out_dir = os.path.join(HEDGEHOG_DIR),
        delimiter = DELIMITER,
        filtering_method = FILTERING_METHOD
    container:
        PYTHON_CONTAINER
    script:
        os.path.join(SRC, "filter_panda.py")

## --------------- ##
## Running BiHiDeF ##
## --------------- ##

rule run_bihidef:
    """
    This rule runs the BiHiDeF algorithm.

    BiHiDeF is available at
    """
    input:
        net = PANDA_NET_FILTERED
    output:
        gene_communities = GENE_COMMUNITIES
    params:
        run_script = os.path.join(SRC, "run_bihidef.py"), \
        max_communities = MAX_COMMUNITIES, \
        max_resolution = MAX_RESOLUTION, \
        output_prefix_reg = REG_TAG, \
        output_prefix_tar = TAR_TAG, \
        outdir = BIHIDEF_RUN_DIR, \
        resource_log = os.path.join(BIHIDEF_RUN_DIR, "run_resources.log"), \
        filtering_method = FILTERING_BIHIDEF,
        ncores = NCORES
    container:
        PYTHON_CONTAINER
    script:
        params.run_script