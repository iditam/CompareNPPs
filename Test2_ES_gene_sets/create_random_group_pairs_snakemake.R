#!/usr/bin/r

# An envelop to call "create_random_group_pairs" from a snakemake file, using the "snakemake" object
source("create_random_group_pairs.R")
create_random_group_pairs(gene_set_size = as.numeric(snakemake@wildcards[["geneSetSize"]]),
                          all_genes_file = snakemake@input[["all_genes_file"]],
                          output_pairs_file = snakemake@output[[1]],
                          config = snakemake@config)
