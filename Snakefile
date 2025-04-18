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
SRC = config["src_dir"]

# Other params
DELIMITER = config["delimiter"]
TISSUE = config["tissue"]

## ----------------- ##
## SiSaNA parameters ##
## ----------------- ##

# SiSaNA directories
SISANA_DIR = os.path.join(OUTPUT_DIR, "{tissue_type}", "sisana")

# SiSaNA inputs
EXP_FILE = os.path.join(DATA_DIR, "{tissue_type}", "gtex_" + "{tissue_type}" + "_exp.tsv")
SISANA_CONFIG = os.path.join(SISANA_DIR, "sisana_params.yml")

# SiSana params
EXP = config["exp_file"]
MOTIF_PRIOR_TAG = config["motif_prior_file"]
PPI_PRIOR_TAG = config["ppi_prior_file"]
MOTIF_PRIOR = os.path.join(DATA_DIR, MOTIF_PRIOR_TAG + ".tsv")
PPI_PRIOR = os.path.join(DATA_DIR, PPI_PRIOR_TAG + ".tsv")

# SiSaNA outputs
EXPRESSION_FILTERED = os.path.join(SISANA_DIR, "preprocess", EXP + "_filtered.txt")
MOTIF_PRIOR_FILTERED = os.path.join(SISANA_DIR, "preprocess", MOTIF_PRIOR_TAG + "_filtered.txt")
PPI_PRIOR_FILTERED = os.path.join(SISANA_DIR, "preprocess", PPI_PRIOR_TAG + "_filtered.txt")
STATS = os.path.join(SISANA_DIR, "preprocess", EXP + "_filtering_statistics.txt")
PANDA_NET = os.path.join(SISANA_DIR, "network", "panda_network.txt")
MOTIF_PRIOR_FILTERED = os.path.join(SISANA_DIR, "preprocess", MOTIF_PRIOR_TAG + "_filtered.txt")

## ---------------- ##
## ELAND parameters ##
## ---------------- ##

# ELAND directories
ELAND_DIR = os.path.join(OUTPUT_DIR, "{tissue_type}", "eland")

# filtering params
PRIOR_ONLY = config["prior_only"]

# network processing outputs
PANDA_EDGELIST = os.path.join(ELAND_DIR, "panda_network_edgelist.txt")
PANDA_NET_FILTERED = os.path.join(ELAND_DIR, "panda_network_filtered.txt")

## ---------------- ------ ##
## Benchmarking parameters ##
## ------------- --------- ##
# directories
BENCHMARK_DIR = os.path.join(ELAND_DIR, "benchmarking")

# params
BENCH_RESOLUTION = config["bench_resolution"]

# outputs
FILTERING_BENCH = os.path.join(BENCHMARK_DIR, "filtering_benchmark_R" + "{bench_resolution}" + ".txt")
FILTERING_BENCH_DF = expand(os.path.join(BENCHMARK_DIR, "filtering_benchmark_df.txt"), tissue_type = TISSUE)
FILTERING_BENCH_PLOT = os.path.join(BENCHMARK_DIR, "filtering_benchmark_plot.pdf")

## ------------------ ##
## BiHiDeF parameters ##
## ------------------ ##

# BiHiDeF params
TAR_TAG = config["target_tag"]
REG_TAG = config["regulator_tag"]
MAX_COMMUNITIES = config["max_communities"]
MAX_RESOLUTION = config["max_resolution"]
MAX_GENES = config["max_genes"]
MIN_GENES = config["min_genes"]

# BiHiDeF directories
BIHIDEF_DIR = os.path.join(ELAND_DIR, "bihidef")
BIHIDEF_RUN_DIR = os.path.join(BIHIDEF_DIR, "C" + str(MAX_COMMUNITIES) + "_R" + str(MAX_RESOLUTION)) # for running with different params

# BiHiDeF outputs
GENE_COMMUNITIES = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG  + ".nodes")
SELECTED_COMMUNITIES = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + "_selected_communities.gmt")
COMMUNITY_STATS = os.path.join(BIHIDEF_RUN_DIR, TAR_TAG + "_community_stats.txt")

## ----------------- ##
## sambar parameters ##
## ----------------- ##

# SAMBAR directories
SAMBAR_DIR = os.path.join(ELAND_DIR, "sambar")
SAMBAR_RUN_DIR = os.path.join(SAMBAR_DIR, "C" + str(MAX_COMMUNITIES) + "_R" + str(MAX_RESOLUTION))

# SAMBAR input files
MUT_DATA = os.path.join(DATA_DIR, config["mut_data"])
ESIZE = os.path.join(DATA_DIR, config["esize_file"])
CAN_GENES = os.path.join(DATA_DIR, config["can_genes"])

# sambar outputs
PATHWAY_SCORES = os.path.join(SAMBAR_RUN_DIR, "pt_out.csv")
MUTATION_SCORES = os.path.join(SAMBAR_RUN_DIR, "mt_out.csv")
DIST_MATRIX = os.path.join(SAMBAR_RUN_DIR, "dist_matrix.csv")

## ------------------------------ ##
## Downstream analysis parameters ##
## ------------------------------ ##

# Analysis directories
ANALYSIS_DIR = os.path.join(OUTPUT_DIR, "{tissue_type}", "analysis")
ANALYSIS_RUN_DIR = os.path.join(ANALYSIS_DIR, "C" + str(MAX_COMMUNITIES) + "_R" + str(MAX_RESOLUTION))

## GO enrichment params ##

# inputs
GENE_BACKGROUND = os.path.join(DATA_DIR, config["bg_file"])
AUTO_BG = config["auto_bg"]
SAVE_ALL = config["save_all"]
SIG_THRESH = config["sig_thresh"]
STATISTIC = config["statistic"]
ALG = config["algorithm"]

# directories
GO_DIR = os.path.join(ANALYSIS_RUN_DIR, "GO_results_" + ALG + "_" + STATISTIC)

# outputs
GO_ENRICHMENT = os.path.join(GO_DIR, STATISTIC + "_GO_summary.txt")

## Clustering params ##

# clustering params
BINARISE = config["binarise"]
LOG_TRANSFORM = config["log_transform"]

if BINARISE:
    TAG = "bin"
elif LOG_TRANSFORM:
    TAG = "log"
else:
    TAG = "raw"

# outputs
CLUST_HEATMAP = os.path.join(ANALYSIS_RUN_DIR, "community_scores_clusters_" + TAG + ".pdf")

## Top mutated pathways ##
# params
N_TOP_COMM = config["n_top_comm"]

# output
TOP_MUTATED_COMMUNITIES = os.path.join(ANALYSIS_RUN_DIR, "top_mut_comm_" + str(N_TOP_COMM) + ".RData")
GENE_MUTATION_SUMMARY = os.path.join(ANALYSIS_RUN_DIR, "gene_mutation_summary_" + str(N_TOP_COMM) + ".txt")
TOP_GENE_SUMMARY = os.path.join(ANALYSIS_RUN_DIR, "top_gene_summary_" + str(N_TOP_COMM) + ".RData")

##-------##
## Rules ##
##-------##

## Rule ALL ##
rule all:
    input:
        expand(PANDA_NET_FILTERED, tissue_type = TISSUE), \
        #expand(GO_ENRICHMENT, tissue_type = TISSUE), \
        #CLUST_HEATMAP, \
        #TOP_MUTATED_COMMUNITIES, \
        #GENE_MUTATION_SUMMARY, \
        #TOP_GENE_SUMMARY, \
        expand(FILTERING_BENCH, tissue_type = TISSUE, bench_resolution = BENCH_RESOLUTION), \
        expand(FILTERING_BENCH_PLOT, tissue_type = TISSUE)

##-----------------------##
## Making PANDA networks ##
##-----------------------##

rule generate_sisana_params:
    """
    This rule generates the SiSaNA config file.
    """
    output:
        config_file = SISANA_CONFIG
    params:
        script = os.path.join(SRC, "utils/generate_sisana_config.py"),
        exp = EXP_FILE,
        motif = MOTIF_PRIOR,
        ppi = PPI_PRIOR,
        number = 5,
        outdir = os.path.join(OUTPUT_DIR, "{tissue_type}", "sisana", "preprocess"),
        processed_paths = os.path.join(OUTPUT_DIR, "{tissue_type}", "sisana", "tmp", "processed_data_paths.yml"),
        method = "panda",
        pandafilepath = os.path.join(OUTPUT_DIR, "{tissue_type}", "sisana", "network", "panda_network.txt")
    container:
        PYTHON_CONTAINER
    message:
        "; Generating SiSaNA config file with script {params.script}"
    shell:
        """
        python {params.script} --exp {params.exp} --motif {params.motif} --ppi {params.ppi} --number {params.number} --outdir {params.outdir} --processed_paths {params.processed_paths} --method {params.method} --pandafilepath {params.pandafilepath} --output {output.config_file}
        """

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
        sisana_config = SISANA_CONFIG
    output:
        EXPRESSION_FILTERED, \
        MOTIF_PRIOR_FILTERED, \
        PPI_PRIOR_FILTERED, \
        STATS, \
        PANDA_NET
    container:
        PYTHON_CONTAINER
    message: 
        "; Running sisana preprocess on {input.sisana_config}."
    shell:
        """
        echo sisana preprocess {input.sisana_config}
        sisana preprocess {input.sisana_config}

        echo sisana generate {input.sisana_config}
        sisana generate {input.sisana_config}
        """

# ##-------------------------------------##
# ## Filtering PANDA network for BiHiDef ##
# ##-------------------------------------##

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
        filtered_net = PANDA_NET_FILTERED, \
        updated_prior_sep = MOTIF_PRIOR_FILTERED.replace(".txt", "_updated_sep.txt"), \
        updated_panda_sep = PANDA_NET.replace(".txt", "_updated_sep.txt") 
    params:
        process = os.path.join(SRC, "process_networks/process_panda.py"), \
        filter = os.path.join(SRC, "process_networks/filter_panda.py"), \
        delimiter = DELIMITER, \
        prior_only = PRIOR_ONLY
    container:
        PYTHON_CONTAINER
    message: 
        "; Processing and filtering PANDA network." \
        "Running {params.process} on {input.panda} and {input.prior} to create {output.edgelist} with --delimiter {params.delimiter}." \
        "Running {params.filter} on {output.updated_panda_sep} and {input.prior} to create {output.filtered_net} with --delimiter {params.delimiter} and --prior_only {params.prior_only}."
    shell:
        """
        python {params.process} {input.panda} {input.prior} {output.edgelist} --delimiter '{params.delimiter}'
        python {params.filter} {output.updated_prior_sep} {output.edgelist} {output.filtered_net} --delimiter '{params.delimiter}' --prior_only '{params.prior_only}'
        """

## ------------------- ##
## Benchmark filtering ##
## ------------------- ##
rule benchmark_filtering:
    """
    This rule benchmarks filtering methods for the PANDA network.

    Inputs
    ------
    PANDA_NET_FILTERED:
        A TXT file with the filtered PANDA network.
    MOTIF_PRIOR_FILTERED:
        A TXT file with the filtered motif prior.
    PANDA_EDGELIST:
        A TXT file with the PANDA edgelist.
    ------
    Outputs
    -------
    BENCHMARK_FILTERED:
        A TXT file with the benchmark data filtered.
    """
    input:
        panda_network_filtered = PANDA_NET_FILTERED, \
        updated_prior_sep = MOTIF_PRIOR_FILTERED.replace(".txt", "_updated_sep.txt"), \
        panda_edgelist = PANDA_EDGELIST
    output:
        filtering_bench = FILTERING_BENCH
    params:
        script = os.path.join(SRC, "process_networks/filter_benchmark.py"), \
        out_dir = os.path.join(BENCHMARK_DIR), \
        delimiter = DELIMITER, \
        resolution = '{bench_resolution}', \
        max_communities = MAX_COMMUNITIES, \
        prior_only = PRIOR_ONLY
    container:
        PYTHON_CONTAINER
    message:
        "; Filtering benchmark data with script {params.script}" \
            "--filtered_net {input.panda_network_filtered}" \
            "--prior_file {input.updated_prior_sep}" \
            "--panda_edgelist {input.panda_edgelist}" \
            "--output_file {output.filtering_bench}" \
            "--delimiter {params.delimiter}" \
            "--resolution {params.resolution}" \
            "--max_communities {params.max_communities}" \
            "--prior_only {params.prior_only}"
    shell:
        """
        python {params.script} --filtered_net {input.panda_network_filtered} --prior_file {input.updated_prior_sep} --panda_edgelist {input.panda_edgelist} --output_file {output.filtering_bench} --delimiter {params.delimiter} --resolution {params.resolution} --max_communities {params.max_communities} --prior_only {params.prior_only}
        """

rule plot_benchmark:
    """
    This rule plots the benchmark data.

    Inputs
    ------
    FILTERING_BENCH_LIST:
        A list of TXT files with the benchmark data for all resolutions.
    ------
    Outputs
    -------
    BENCHMARK_PLOT:
        A PDF file with the benchmark plot.
    """
    input:
        filtering_bench_list=expand(FILTERING_BENCH, tissue_type="{tissue_type}", bench_resolution=BENCH_RESOLUTION)
    output:
        benchmark_plot=FILTERING_BENCH_PLOT
    params:
        script=os.path.join(SRC, "analysis", "plot_filtering_bench.R"), \
        out_dir=os.path.join(BENCHMARK_DIR), \
        tissue_type="{tissue_type}", \
        data_frame=os.path.join(BENCHMARK_DIR, "filtering_benchmark_df.txt"), \
        plot_type="all", \
        plot_title="Filtering benchmark for {tissue_type}", \
        files=lambda wildcards, input: ",".join(input.filtering_bench_list)  # Preprocess the input list
    container:
        ANALYSIS_CONTAINER
    message:
        "; Plotting benchmark data with script {params.script} for tissue {wildcards.tissue_type}"
    shell:
        """
        Rscript {params.script} \
            --files "{params.files}" \
            --output-dir {params.out_dir} \
            --tissue-type {params.tissue_type} \
            --data-frame {params.data_frame} \
            --plot-type {params.plot_type} \
            --plot-title "{params.plot_title}" \
            --plot-file {output.benchmark_plot}
        """

# --------------- ##
# Running BiHiDeF ##
# --------------- ##

rule run_bihidef:
    """
    This rule runs the BiHiDeF algorithm.

    BiHiDeF is available at
    """
    input:
        PANDA_NET_FILTERED
    output:
        gene_communities = GENE_COMMUNITIES
    params:
        run_script = os.path.join(SRC, "eland/run_bihidef.py"), \
        measure_script = os.path.join(SRC, "utils/measure_resources.py"), \
        max_communities = MAX_COMMUNITIES, \
        max_resolution = MAX_RESOLUTION, \
        output_prefix_reg = REG_TAG, \
        output_prefix_target = TAR_TAG, \
        out_dir = BIHIDEF_RUN_DIR, \
        log_file = os.path.join(BIHIDEF_RUN_DIR, "run_log.log")
    container:
        PYTHON_CONTAINER
    message:
        "; Running BiHiDeF on {input} with params:" \
            "--comm_mult {params.max_communities}" \
            "--max_res {params.max_resolution}" \
            "--output_dir {params.out_dir}" \
            "--output_prefix_reg {params.output_prefix_reg}" \
            "--output_prefix_tar {params.output_prefix_target}"
    shell:
        """
        mkdir -p {params.out_dir}
        python {params.measure_script} {params.log_file} "python {params.run_script} {input} --comm_mult {params.max_communities} --max_res {params.max_resolution} \
        --output_dir {params.out_dir} --output_prefix_reg {params.output_prefix_reg} --output_prefix_tar {params.output_prefix_target}"
        """

# # ## --------------------- ##
# # ## Selecting communities ##
# # ## --------------------- ##

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
        script = os.path.join(SRC, "eland/select_communities.py"), \
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

# ## -------------- ##
# ## Running sambar ##
# ## -------------- ##

# rule run_sambar:
#     """
#     This rule runs the sambar algorithm.

#     Sambar is available at:
#     """
#     input:
#         gmt = SELECTED_COMMUNITIES, \
#         mut = MUT_DATA, \
#         esize = ESIZE, \
#         can_genes = CAN_GENES
#     output:
#         pathway_scores = PATHWAY_SCORES, \
#         mutation_scores = MUTATION_SCORES, \
#         dist_matrix = DIST_MATRIX

#     params:
#         script = os.path.join(SRC, "eland/run_sambar.py"), \
#         out_dir = SAMBAR_RUN_DIR
#     container:
#         PYTHON_CONTAINER
#     message:
#         "; Running sambar on {input} with params:" \
#             "--output-dir {params.out_dir} " \
#             "--gmt-file {input.gmt} " \
#             "--mut-file {input.mut} " \
#             "--esize-file {input.esize} " \
#             "--can-genes {input.can_genes} "
#     shell:
#         """
#         mkdir -p {params.out_dir}
#         python {params.script} --gmt-file {input.gmt} --mut-file {input.mut} --esize-file {input.esize} --output-dir {params.out_dir} --can-genes {input.can_genes}
#         """

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
        gmt = SELECTED_COMMUNITIES, \
        bg = GENE_BACKGROUND
    output:
        go_enrichment = GO_ENRICHMENT
    params:
        script = os.path.join(SRC, "analysis/GO_enrichment.R"), \
        auto_bg = AUTO_BG, \
        out_dir = GO_DIR, \
        save_all = SAVE_ALL, \
        sig_thresh = SIG_THRESH, \
        statistic = STATISTIC, \
        algorithm = ALG

    container:
        ANALYSIS_CONTAINER
    message:
        "; Running GO enrichment with script {params.script}" \
            "--gmt-file {input.gmt} " \
            "--bg-file {input.bg} " \
            "--out-dir {params.out_dir} " \
            "--auto-bg {params.auto_bg} " \
            "--save-all {params.save_all} " \
            "--thresh {params.sig_thresh} " \
            "--statistic {params.statistic} " \
            "--algorithm {params.algorithm} "
    shell:
        """
        mkdir -p {params.out_dir}
        echo Rscript {params.script} --gmt-file {input.gmt} --bg-file {input.bg} --auto-bg {params.auto_bg} --save-all {params.save_all} --sig-thresh {params.sig_thresh} --statistic {params.statistic} --algorithm params.algorithm --output-dir {params.out_dir}
        Rscript {params.script} --gmt-file {input.gmt} --bg-file {input.bg} --auto-bg {params.auto_bg} --save-all {params.save_all} --thresh {params.sig_thresh} --statistic {params.statistic} --algorithm {params.algorithm} --output-dir {params.out_dir}

        """

# ## ----------------------- ##
# ## Top mutated communities ##
# ## ----------------------- ##

# rule top_mutated_communities:
#     """
#     This rule extracts the top n most mutated communities and finds the top genes in these communities that are mutated.

#     Input:
#     ------
#     SELECTED_COMMUNITIES:
#         A GMT file with the selected communities.
#     MUT_SCORES:
#         File with mutation scores.
#     PATHWAY_SCORES:
#         File with pathways scores.
    
#     Output:
#     -------
#     TOP_MUTATED_COMMUNITIES:
#         RData file with the top mutated communities.
#     GENE_MUTATION_SUMMARY:
#         A TXT file with the gene mutation summary.
#     TOP_GENE_SUMMARY:
#         RData file with the top gene summary.
#     """
#     input:
#         communities = SELECTED_COMMUNITIES, \
#         mut = MUT_DATA, \
#         pathway_scores = PATHWAY_SCORES
#     output:
#         top_mutated_communities = TOP_MUTATED_COMMUNITIES, \
#         gene_mutation_summary = GENE_MUTATION_SUMMARY, \
#         top_gene_summary = TOP_GENE_SUMMARY
#     container:
#         ANALYSIS_CONTAINER
#     params:
#         script = os.path.join(SRC, "analysis/top_mut_comm.R"), \
#         out_dir = ANALYSIS_RUN_DIR, \
#         n_top_comm = N_TOP_COMM
#     message:
#         """
#         ; Running Rscript {params.script} --communities {input.communities} --mut {input.mut} --scores {input.pathway_scores} --output-dir {params.out_dir} --n-comm {params.n_top_comm}
#         """
#     shell:
#         """
#         mkdir -p {params.out_dir}
#         Rscript {params.script} --communities {input.communities} --mut {input.mut} --scores {input.pathway_scores} --output-dir {params.out_dir} --n-comm {params.n_top_comm}
#         """

# ## --------------------------- ##
# ## Clustering community scores ##
# ## --------------------------- ##

# rule cluster_scores:
#     """
#     This rule clusters the communitiy scores obtained from SAMBAR.

#     Inputs
#     ------
#     PATHWAY_SCORES:
#         csv file with the community scores as output by sambar.
#     DIST_MATRIX:
#         csv file with the distance matrix as output by sambar. 

#     Outputs
#     -------
#     HIERARCHICAL_CLUSTERING:
#         Hierarchical clustering of community scores.
#     CLUST_HEATMAP:
#         A heatmap of the clustering. 
#     """
#     input:
#         community_scores = PATHWAY_SCORES, \
#         dist_matrix = DIST_MATRIX
#     output:
#         heatmap = CLUST_HEATMAP
#     container:
#         ANALYSIS_CONTAINER
#     params:
#         script = os.path.join(SRC, "analysis/cluster_scores.R"), \
#         out_dir = ANALYSIS_RUN_DIR, \
#         binarise = BINARISE, \
#         log_transform = LOG_TRANSFORM
#     message:
#         """
#         ; Running Rscript {params.script} --file {input.community_scores} --output-dir {params.out_dir} --dist {input.dist_matrix} --binarise {params.binarise} --log-transform {params.log_transform}
#         """
#     shell:
#         """
#         echo Rscript {params.script} --file {input.community_scores} --output-dir {params.out_dir} --dist {input.dist_matrix} --binarise {params.binarise} --log-transform {params.log_transform}
#         mkdir -p {params.out_dir}
#         Rscript {params.script} --file {input.community_scores} --output-dir {params.out_dir} --dist {input.dist_matrix} --binarise {params.binarise} --log-transform {params.log_transform}
#         """