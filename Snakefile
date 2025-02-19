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