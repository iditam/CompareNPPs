#!/usr/bin/r

# An envelop to call "calculate_ES" from a snakemake file, using the "snakemake" object
source("calculate_ES.R")
# input_correlation_file, input_gene_group_pairs, output_ES_file, nppI
calculate_ES(input_correlation_file = snakemake@input[["correlation_file"]], 
             input_gene_group_pairs = snakemake@input[["gene_group_pairs"]],
             output_ES_file = snakemake@output[[1]], 
             nppID = snakemake@wildcards[["NPP"]],
             config = snakemake@config)
