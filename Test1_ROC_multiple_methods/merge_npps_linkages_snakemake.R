#!/usr/bin/r

# An envelop to call "merge_npps_linkages" from a snakemake file, using the "snakemake" object
source("merge_npps_linkages.R")
merge_npps_linkages(linkage_file = snakemake@input[["linkage_file"]], 
                    npp_linkage_files = snakemake@input[["npp_linkage_files"]],
                    output_file = snakemake@output[["filename"]],
                    config = snakemake@config)
