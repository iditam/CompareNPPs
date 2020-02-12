#!/usr/bin/r

# An envelop to call "create_random_group_pairs" from a snakemake file, using the "snakemake" object
source("calculate_ES_for_random_groups_partial.R")
calculate_ES_for_random_groups_partial(correlation_file = snakemake@input[["correlation_file"]],
                                       random_pairs_file = snakemake@input[["random_pairs_file"]],
                                       current_range = snakemake@params[["current_range"]],
                                       output_ES_file = snakemake@output[[1]])
