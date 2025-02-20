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

# Directories
DATA_dir = config["data_dir"]
SISANA_DIR = config["sisana_dir"]
SISANA_OUTPUT_DIR = os.path.join(SISANA_DIR, "output") # might want to change later as user can change output dir in the sisana params, but it'll do for now

# Other params
EXP = config["exp_file"]
MOTIF_PRIOR = config["motif_prior_file"]
PPI_PRIOR = config["ppi_prior_file"]

# Input files
SISANA_CONFIG = os.path.join(SISANA_DIR, config["sisana_config"])

# Output files
EXPRESSION_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", EXP + "_filtered.txt")
MOTIF_PRIOR_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", MOTIF_PRIOR + "_filtered.txt")
PPI_PRIOR_FILTERED = os.path.join(SISANA_OUTPUT_DIR, "preprocess", PPI_PRIOR + "_filtered.txt")
STATS = os.path.join(SISANA_OUTPUT_DIR, "preprocess", EXP + "_filtering_statistics.txt")
PANDA_NET = os.path.join(SISANA_OUTPUT_DIR, "network", "panda_network.txt")


## Rule ALL##
rule all:
    input: 
        EXPRESSION_FILTERED, \
        MOTIF_PRIOR_FILTERED, \
        PPI_PRIOR_FILTERED, \
        STATS, \
        PANDA_NET

## Rules ##
rule run_sisana:
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
    PANDA_NET:
        A TXT file with the PANDA network.
    """
    input:
        SISANA_CONFIG
    output:
        EXPRESSION_FILTERED, \
        MOTIF_PRIOR_FILTERED, \
        PPI_PRIOR_FILTERED, \
        STATS, \
        PANDA_NET

    message: 
        "; Running sisana preprocess on {input}."
    shell:
        """
        sisana preprocess {input}
        sisana generate {input}
        """
