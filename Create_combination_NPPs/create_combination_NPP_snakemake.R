#!/usr/bin/r

# An envelop to call "create_combination_NPP" from a snakemake file, using the "snakemake" object
source("create_combination_NPP.R")
create_combination_NPP(threshold =  snakemake@wildcards[["threshold"]],
                       transformation =  snakemake@wildcards[["transformation"]],
                       divide_by_query =  snakemake@wildcards[["divide_by_query"]],
                       colzscore = snakemake@wildcards[["colzscore"]],
                       npp_filename = snakemake@output[["npp_filename"]])
