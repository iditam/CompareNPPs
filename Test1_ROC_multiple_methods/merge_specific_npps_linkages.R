#!/usr/bin/r
library(data.table)
library(stringr)
library(stringi)
library(plyr)

correlation_column = "correlation"

linkage_file <- snakemake@input[["linkage_file"]]
npp_linkage_files <- snakemake@input[["npp_linkage_files"]]
output_file <- snakemake@output[["filename"]]
config <- snakemake@config

linkages <- fread(input = linkage_file, verbose = F, showProgress = F)
npp_linkages <- lapply(npp_linkage_files, function(fname) { fread(input = fname, verbose = F, showProgress = F) })

# set the names(npp_linkages) as the identifications of the NPP versions being compared
# will use these names as the column titles in the merged linkages file
names(npp_linkages) <- sub(".*\\-(.+\\-.+).csv",'\\1', npp_linkage_files) 

for(npp_name in names(npp_linkages)) {
  setnames(npp_linkages[[npp_name]], correlation_column, npp_name)
  linkages <- merge(linkages, npp_linkages[[npp_name]], by = c("LinkID","gene_pair"), all.x = T)
}

fwrite(linkages, file = output_file)

