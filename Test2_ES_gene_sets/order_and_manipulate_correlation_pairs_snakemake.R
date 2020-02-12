#!/usr/bin/r

# An envelop to call "order_and_zero_correlation_pairs" from a snakemake file, using the "snakemake" object
source("order_and_manipulate_correlation_pairs.R")
order_and_manipulate_correlation_pairs(correlation_file = snakemake@input[["correlation_file"]],
                                       output_correlation_file = snakemake@output[[1]],
                                       config = snakemake@config)
