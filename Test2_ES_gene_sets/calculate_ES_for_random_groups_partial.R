library(data.table)
# library(parallel)

# OUTPUT:
# a file with ES values of the given nppID, enrichment for all the randomized groups of the given group size


# the Enrichment Score function, based on code from GSEA ----
GSEA.EnrichmentScore <- function(gene_set_pairs, correlation_pairs_NPP) {
  N <- nrow(correlation_pairs_NPP)
  Nh <- length(gene_set_pairs)
  Nm <-  N - Nh
  tag.indicator <- sign(match(correlation_pairs_NPP$gene_pair, gene_set_pairs, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  correl.vector <- correlation_pairs_NPP$correlation
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  return(signif(max.ES, digits = 3))
  
  # correlation_pairs_NPP[, gene_pair_found := 0]
  # correlation_pairs_NPP[gene_pair %in% gene_set_pairs, gene_pair_found := 1]
  # 
  # if (weighted.score.type == 0) {
  #   correlation_pairs_NPP[ ,correl.vector:=1]
  #   sum.correl.tag <- Nh
  # }
  # else if (weighted.score.type == 1) {
  #   correlation_pairs_NPP[ ,correl.vector:=0] # from former call to this function
  #   # add correl.vector with values only for gene_pair that appear in gene_set_pairs.
  #   # other rows will have correl.vector=NA
  #   correlation_pairs_NPP[gene_pair_found==1, correl.vector:=abs(correlation)]
  #   sum.correl.tag <- correlation_pairs_NPP[, sum(correl.vector)]
  # }
  # else {
  #   correlation_pairs_NPP[ ,correl.vector:=0]  # from former call to this function
  #   # add correl.vector with values only for gene_pair that appear in gene_set_pairs.
  #   # other rows will have correl.vector==0
  #   correlation_pairs_NPP[gene_pair_found==1, correl.vector:=abs(correlation^weighted.score.type)]
  #   sum.correl.tag <- correlation_pairs_NPP[, sum(correl.vector)]
  # }
  # norm.tag    <- 1.0/sum.correl.tag
  # minus.norm.no.tag <- -1.0/Nm # changed to minus!
  # correlation_pairs_NPP[ ,value_for_cumsum:=minus.norm.no.tag][gene_pair_found==1,value_for_cumsum:=correl.vector*norm.tag]
  # correlation_pairs_NPP[ ,RES := cumsum(value_for_cumsum)]
  # max.ES <- correlation_pairs_NPP[ ,max(RES)]
  # # changed to take only max.ES
  # # min.ES <- correlation_pairs_NPP[ ,min(RES)]
  # # if (max.ES > - min.ES) {
  # #   ES <- signif(max.ES, digits = 5)
  # # } else {
  # #   ES <- signif(min.ES, digits=5)
  # # }
  # # return(ES)
  # return(signif(max.ES, digits = 3))
}

calculate_ES_for_random_groups_partial <- function(correlation_file, random_pairs_file, 
                                                   current_range, output_ES_file) {
 
  randomPairs <- fread(random_pairs_file, header = F, verbose = F, showProgress = F)

  # load correlation pairs
  correlation_pairs_NPP <- fread(input = correlation_file, verbose = F, showProgress = F)


 
  # each row in randomPairs is actually a vector of all gene pairs of the random gene set
  ES_values_vec <- sapply(seq(current_range[1], current_range[2]),
                             function(i) {
                               GSEA.EnrichmentScore(gene_set_pairs = randomPairs[i],
                                                    correlation_pairs_NPP)
                             })

  writeLines(as.character(ES_values_vec), con = output_ES_file)
}
