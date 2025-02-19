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

# Input files
SISANA_CONFIG = config["sisana_config"]

# Output files
EXPRESSION_FILTERED = os.path.join(SISANA_DIR, "/output/preprocess/", "{filename}_filtered.txt")
MOTIF_PRIOR_FILTERED = os.path.join(SISANA_DIR, "/output/preprocess/", "{filename}_filtered.txt")
PPI_PRIOR_FILTERED = os.path.join(SISANA_DIR, "/output/preprocess/", "{filename}_filtered.txt")
STATS = os.path.join(SISANA_DIR, "/output/preprocess/", "{filename}_filtering_statistics.txt")

# Other params
EXP = config["exp_file"]
MOTIF_PRIOR = config["motif_prior_file"]
PPI_PRIOR = config["ppi_prior_file"]

## Rule ALL##
rule all:
    input: 
        expand(EXPRESSION_FILTERED, filename = EXP), \
        expand(MOTIF_PRIOR_FILTERED, filename = MOTIF_PRIOR), \
        expand(PPI_PRIOR_FILTERED, filename = PPI_PRIOR), \
        expand(STATS, filename = EXP)

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
    """
    input:
        sisana_config = SISANA_CONFIG
    output:
        EXPRESSION_FILTERED, \
        MOTIF_PRIOR_FILTERED, \
        PPI_PRIOR
    message: 
        "; Running sisana preprocess on {input}."
    params:

    shell:
        """
        sisana preprocess {input.sisana_config}
        sisana generate {input.sisana_config}
        """
