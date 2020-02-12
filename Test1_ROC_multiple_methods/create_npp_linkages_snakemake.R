#!/usr/bin/r

# An envelop to call "create_npp_linkages" from a snakemake file, using the "snakemake" object
source("create_npp_linkages.R")
create_npp_linkages(linkage_file = snakemake@input[["linkage_file"]], 
                    correlation_file = snakemake@input[["correlation_file"]],
                    output_file = snakemake@output[["filename"]])
