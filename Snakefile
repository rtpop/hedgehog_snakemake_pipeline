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


##-------##
## Rules ##
##-------##

## Rule ALL ##
rule all:
    input:
        expand(PANDA_NET_FILTERED, tissue_type = TISSUE)

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
    MOTIF_PRIOR:
        A TXT file with filtered motif prior.
    ------
    Outputs
    -------
    # PANDA_EDGELIST:
    #     A TXT file with the PANDA network processed as an edgelist compatible with networkx.
    PANDA_NET_FILTERED:
        A TXT file with the PANDA edgelist filtered.
    """
    input:
        panda = os.path.join(DATA_DIR, "{tissue_type}", "panda_network_edgelist.txt"), \
        prior = MOTIF_PRIOR
    output:
        filtered_net = PANDA_NET_FILTERED
    params:
        filter = os.path.join(SRC, "process_networks/filter_panda.py"), \
        delimiter = DELIMITER, \
        prior_only = PRIOR_ONLY, \
        out_dir = ELAND_DIR
    container:
        PYTHON_CONTAINER
    message: 
        "; Processing and filtering PANDA network." \
        # "Running {params.process} on {input.panda} and {input.prior} to create {output.edgelist} with --delimiter {params.delimiter}." \
        "Running {params.filter} on {input.panda} and {input.prior} to create {output.filtered_net} with --delimiter {params.delimiter} and --prior_only {params.prior_only}."
    shell:
        """
        mkdir -p {params.out_dir}
        python {params.filter} {input.prior} {input.panda} {output.filtered_net} --delimiter '{params.delimiter}' --prior_only '{params.prior_only}'
        """

# ## ------------------- ##
# ## Benchmark filtering ##
# ## ------------------- ##
rule benchmark_filtering:
    """
    This rule benchmarks filtering methods for the PANDA network.

    Inputs
    ------
    PANDA_NET_FILTERED:
        A TXT file with the filtered PANDA network.
    MOTIF_PRIOR:
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
        prior = MOTIF_PRIOR, \
        panda_edgelist = os.path.join(DATA_DIR, "{tissue_type}", "panda_network_edgelist.txt")
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
            "--prior_file {input.prior}" \
            "--panda_edgelist {input.panda_edgelist}" \
            "--output_file {output.filtering_bench}" \
            "--delimiter {params.delimiter}" \
            "--resolution {params.resolution}" \
            "--max_communities {params.max_communities}" \
            "--prior_only {params.prior_only}"
    shell:
        """
        python {params.script} --filtered_net {input.panda_network_filtered} --prior_file {input.prior} --panda_edgelist {input.panda_edgelist} --output_file {output.filtering_bench} --delimiter {params.delimiter} --resolution {params.resolution} --max_communities {params.max_communities} --prior_only {params.prior_only}
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
        files=lambda wildcards, input: ",".join(input.filtering_bench_list), \
        include_unfiltered=UNFILTERED
    container:
        ANALYSIS_CONTAINER
    message:
        "; Plotting benchmark data with script {params.script} for tissue {wildcards.tissue_type} with params:" \
            "--files {params.files}" \
            "--output-dir {params.out_dir}" \
            "--tissue-type {params.tissue_type}" \
            "--data-frame {params.data_frame}" \
            "--plot-type {params.plot_type}" \
            "--plot-title {params.plot_title}"
    shell:
        """
        Rscript {params.script} \
            --files "{params.files}" \
            --output-dir {params.out_dir} \
            --tissue-type {params.tissue_type} \
            --data-frame {params.data_frame} \
            --plot-type {params.plot_type} \
            --plot-title "{params.plot_title}" \
            --plot-file {output.benchmark_plot} \
            --include-unfiltered {params.include_unfiltered}
        """

rule consolidate_benchmark:
    """
    This rule consolidates the benchmark data into a single file.

    Inputs
    ------
    FILTERING_BENCH_LIST:
        A list of TXT files with the benchmark data for all resolutions.
    ------
    Outputs
    -------
    FILTERING_BENCH_DF:
        A TXT file with the consolidated benchmark data.
    """
    input:
        filtering_bench_list=expand(FILTERING_BENCH_DF, tissue_type=TISSUE)
    output:
        filtering_bench_df=FILTERING_BENCH_CONSOLIDATED
    params:
        script=os.path.join(SRC, "utils/consolidate_benchmark.R"), \
        files=lambda wildcards, input: ",".join(input.filtering_bench_list)
    container:
        ANALYSIS_CONTAINER
    message:
        "; Consolidating benchmark data with script {params.script} with params:" \
            "--files {params.files}" \
            "--output-file {output.filtering_bench_df}"
    shell:
        """
        Rscript {params.script} --files "{params.files}" --output {output.filtering_bench_df}
        """

rule plot_bench_heatmap:
    """
    This rule plots the benchmark heatmap.

    Inputs
    ------
    FILTERING_BENCH_CONSOLIDATED:
        A TXT file with the consolidated benchmark data.
    ------
    Outputs
    -------
    BENCHMARK_HEATMAP:
        A PDF file with the benchmark heatmap.
    """
    input:
        filtering_bench_df = FILTERING_BENCH_CONSOLIDATED
    output:
        benchmark_heatmap = FILTERING_HEATMAP
    params:
        script = os.path.join(SRC, "analysis/plot_filtering_bench_heatmap.R"), \
        metric = "Modularity",  \
        filtering_method = "Prior filtered"
    container:
        ANALYSIS_CONTAINER
    message:
        "; Plotting benchmark heatmap with script {params.script}"
    shell:
        """
        Rscript {params.script} --input "{input.filtering_bench_df}" --output {output.benchmark_heatmap} --metric {params.metric} --filtering "{params.filtering_method}"
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

rule consolidate_community_stats:
    """
    This rule consolidates the community stats into a single file.

    Inputs
    ------
    COMMUNITY_STATS:
        A TXT file with the community stats.
    ------
    Outputs
    -------
    COMMUNITY_STATS:
        A TXT file with the consolidated community stats.
    """
    input:
        community_stats = expand(COMMUNITY_STATS, tissue_type = TISSUE), \
        selected_communities = expand(SELECTED_COMMUNITIES, tissue_type = TISSUE), \
        go_summery = expand(GO_ENRICHMENT, tissue_type = TISSUE)
    output:
        consolidated_community_stats = COMMUNITY_STATS_CONSOLIDATED
    params:
        script = os.path.join(SRC, "utils/consolidate_community_stats.py"), \
        stats_joined = lambda wildcards, input: ",".join(input.community_stats), \
        gmt_joined = lambda wildcards, input: ",".join(input.selected_communities), \
        go_summary_joined = lambda wildcards, input: ",".join(input.go_summery)
    container:
        PYTHON_CONTAINER
    message:
        "; Consolidating community stats with script {params.script} " \
            "{input.community_stats} " \
            "{input.selected_communities} " \
            "{output.consolidated_community_stats}"
    shell:
        """
        python {params.script} "{params.stats_joined}" "{params.gmt_joined}" "{params.go_summary_joined}" {output.consolidated_community_stats}
        """

rule plot_n_communities:
    """
    This rule plots the number of selected communities.

    Inputs
    ------
    COMMUNITY_STATS_CONSOLIDATED:
        A TXT file with the consolidated community stats.
    ------
    Outputs
    -------
    COMMUNITY_PLOT:
        A PDF file with the community plot.
    """
    input:
        df = COMMUNITY_STATS_CONSOLIDATED
    output:
        community_plot = COMMUNITY_PLOT
    params:
        script = os.path.join(SRC, "analysis/community_plot.R")
    container:
        ANALYSIS_CONTAINER
    message:
        "; Plotting communities with script {params.script}"
    shell:
        """
        Rscript {params.script} --input {input.df} --output {output.community_plot}
        """

rule plot_n_enriched_communities:
    """
    This rule plots the number of enriched communities.

    Inputs
    ------
    GO_ENRICHMENT:
        A TXT file with the GO enrichment results.
    ------
    Outputs
    -------
    GO_ENRICHMENT_PLOT:
        A PDF file with the GO enrichment plot.
    """
    input:
        go_enrichment = GO_ENRICHMENT
    output:
        go_enrichment_plot = GO_ENRICHMENT_BARPLOT
    params:
        script = os.path.join(SRC, "analysis/go_enrichment_plot.R"), \
        out_dir = GO_DIR
    container:
        ANALYSIS_CONTAINER
    message:
        "; Running GO enrichment plot with script {params.script}" \
            "--go-enrichment {input.go_enrichment} " \
            "--out-dir {params.out_dir} " \
            "--plot-file {output.go_enrichment_plot} "
    shell:
        """
        mkdir -p {params.out_dir}
        Rscript {params.script} --input {input.go_enrichment} --out-dir {params.out_dir} --file-name {output.go_enrichment_plot}
        """