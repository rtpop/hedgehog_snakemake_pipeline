# Run from directory containing Snakefile
# snakemake --cores 1 -np ### For dry run

# Libraries
import os 
import sys
import glob
from pathlib import Path
import time

# Config
global CONFIG_PATH
CONFIG_PATH = "config.yaml"
configfile: CONFIG_PATH

# Containers
SISANA_CONTAINER = config["sisana_container"]
ELAND_CONTAINER = config["eland_container"]
ANALYSIS_CONTAINER = config["analysis_container"]

# Params
SISANA_CONFIG = config["sisana_config"]

# Input files

# Output files

## Rule ALL##
rule all:
    input: 

## Rules ##
rule run_sisana
    """
    This rule runs PANDA using the SiSaNA pipeline.

    SiSaNA is available at
    https://github.com/kuijjerlab/sisana

    Inputs
    ------
    SISANA_CONFIG:
        Config yml file for SiSaNA.
    ------
    Outputs
    -------
    EXPRESSION_FILTERED:
        A TXT file with filtered expression.
    MOTIF_PRIOR_FILTERED:
        A TXT file with filtered motif prior.
    PPI_PRIOR_FILTERED:
        A TXT file with filtered PPI prior.
    STATS:
        A file containing information about genes filtered.    
    """
    input:
        SISANA_CONFIG
    output:
        EXPRESSION_FILTERED, \
        MOTIF_PRIOR_FILTERED, \
        PPI_PRIOR
    message: 
        "; Running sisana preprocess on {input}."
    params:

    shell:
        """
        sisana preprocess {input}
        """
