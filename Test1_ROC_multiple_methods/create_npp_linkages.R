#!/usr/bin/r

library(data.table)

create_npp_linkages <- function(linkage_file, correlation_file, output_file) {
   
  linkages <- fread(input = linkage_file, verbose = F, showProgress = F)
  linkages <- linkages[ ,.(LinkID,gene_pair)]
  correlation_pairs_NPP <- fread(input = correlation_file, verbose = F, showProgress = F)
  
  setkey(correlation_pairs_NPP, gene_pair)
  setkey(linkages, gene_pair)
  npp_linkages <- merge(linkages, correlation_pairs_NPP, all.x=TRUE)
  # npp_linkages[is.na(correlation)]
  
  setkey(npp_linkages, LinkID)
  fwrite(npp_linkages, file = output_file)
}
