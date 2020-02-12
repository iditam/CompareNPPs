#!/usr/bin/r
library(data.table)
library(stringr)
library(stringi)
library(plyr)

correlation_column = "correlation"

merge_npps_linkages <- function(linkage_file, npp_linkage_files, output_file, config) {

  linkages <- fread(input = linkage_file, verbose = F, showProgress = F)
  npp_linkages <- lapply(npp_linkage_files, function(fname) { fread(input = fname, verbose = F, showProgress = F) })

  # set the names(npp_linkages) as the identifications of the NPP versions being compared
  # will use these names as the column titles in the merged linkages file
  if (("CORRELATIONs" %in% names(config)) && (length(config[["CORRELATIONs"]]) > 1)) {
    names(npp_linkages) <- sub(".+\\-(.+\\-.+)\\.csv",'\\1', npp_linkage_files) # will get the npp method and the correlation, e.g. "UniProt_062018__spearman"
  } else {
    names(npp_linkages) <- sub(".*\\-(.+)\\-.*",'\\1', npp_linkage_files) 
  }

  for(npp_name in names(npp_linkages)) {
    setnames(npp_linkages[[npp_name]], correlation_column, npp_name)
    linkages <- merge(linkages, npp_linkages[[npp_name]], by = c("LinkID","gene_pair"), all.x = T)
  }
  
  fwrite(linkages, file = output_file)
  
}
  
