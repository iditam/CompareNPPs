#!/usr/bin/r

# An envelop to call "plot_ROC_NPP_methods" from a snakemake file, using the "snakemake" object
source("ROC_NPP_methods.R")
ROC_NPP_methods(pos_linkages_file = snakemake@input[["positive_linkages"]], 
                neg_linkages_file = snakemake@input[["negative_linkages"]],
                ROC_output_file = snakemake@output[["ROC_filename"]],
                PR_output_file = snakemake@output[["PR_filename"]],
                AUCs_output_file = snakemake@output[["AUCs_filename"]],
                config = snakemake@config)
