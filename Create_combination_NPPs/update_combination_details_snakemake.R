#!/usr/bin/r

# An envelop to call "update_combination_details" from a snakemake file, using the "snakemake" object
source("update_combination_details.R")
update_combination_details(input_list_file = snakemake@input[["filename"]],
                           output_list_file = snakemake@output[["list_file"]],
                           output_details_file = snakemake@output[["details_file"]])
