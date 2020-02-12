#!/usr/bin/r

##########################################################
# Input: 
# input_correlation_file - a file that contains columns:
# correlation,gene_pair
# ordered by correlation descending, and each gene pair is ordered alphabetically ascending (gene1 < gene2)
# config - a list of configuration definitions
#
#  assumes existence of a file geneGroupsPairs_<collection title>.rds
#  that contains all gene set pairs of the given collection.
#  (Created by CompareNPP_3_create_gene_groups_pairs.R)
#
# Action:
#   For each gene set in the collection, calculate its Enrichment Score (ES)
#   in a GSEA-like method, as related to the gene pairs list with correlations (Input 1.)
#
# Output:
#  A file with name <output_ES_file>
#  that contains a table with ES value for each gene set. columns: gene_set,ES
##########################################################
library(data.table)
library(stringr)

source("../CommonCode/ReadGmt.R")


# the Enrichment Score function, based on code from GSEA ----
GSEA.EnrichmentScore <- function(gene_set_name, gene_set_pairs, correlation_pairs_NPP,
                                 nppID, write_RES_to_file, RES_file_name) {  
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
  if (write_RES_to_file)
    writeLines(RES, con = RES_file_name)
  max.ES <- max(RES)
  return(signif(max.ES, digits = 3))
}

# the Enrichment Score function, based on code from GSEA ----
GSEA.EnrichmentScore.per_ranking <- function(gene_set_name, gene_set_pairs, correlation_pairs_NPP,
                                 nppID, write_RES_to_file, RES_file_name, ranking_multiply) {  
  N <- nrow(correlation_pairs_NPP) 
  Nh <- length(gene_set_pairs) 
  Nm <-  N - Nh 
  correlation_pairs_NPP[, gene_pair_found := 0]
  correlation_pairs_NPP[gene_pair %in% gene_set_pairs, gene_pair_found := 1]
  
  sum_rank_for_step <- correlation_pairs_NPP[, sum(rank_for_step)] 
  
  norm.tag    <- 1.0/sum_rank_for_step
  minus.norm.no.tag <- -1.0/Nm # minus because it's a step downward
  correlation_pairs_NPP[ ,value_for_cumsum:=minus.norm.no.tag][gene_pair_found==1,value_for_cumsum:=ranking_multiply*rank_for_step*norm.tag]
  correlation_pairs_NPP[ ,RES := cumsum(value_for_cumsum)]
  
  if (write_RES_to_file)
    fwrite(correlation_pairs_NPP[ ,.(RES)], RES_file_name,
           showProgress = F, verbose = F)
  
  max.ES <- correlation_pairs_NPP[ ,max(RES)]
  return(signif(max.ES, digits = 3))
}

calculate_ES <- function(input_correlation_file, input_gene_group_pairs, 
                         output_ES_file, nppID, config) {
  

  geneGroupsFolder = file.path(config[["geneGroupsRootFolder"]], config[["geneCollectionResource"]])
  collectionName <- config[["collectionName"]]
  outputFolder <- config[["outputFolder"]]
  geneGroups <- get_genes_groups_list(geneGroupsFolder,collectionName)
  
  geneGroupsPairs <- readRDS(input_gene_group_pairs)
  
  # load correlation pairs
  correlation_pairs_NPP <- fread(input = input_correlation_file, verbose = F, showProgress = F)
  
  # order correlation pairs by correlation (descending)
  setorder(correlation_pairs_NPP, -correlation)
  
  if("suffix" %in% names(config))
    if(config[["suffix"]] == "_zeroed)")  {
      # zero negative correlations
      correlation_pairs_NPP[correlation < 0, correlation:=0]
    } else if (str_detect(config[["suffix"]] , "_offset_add")) {
      offset_val <- as.numeric(str_match(config[["suffix"]] , "_offset_add(.+)")[2])
      correlation_pairs_NPP[ , correlation:=correlation+offset_val]
    }
  
  
  
  
  write_RES_to_file <- FALSE
  if(! is.null(config[["write_RES_to_file"]]) && config[["write_RES_to_file"]]) {
    write_RES_to_file <- TRUE
    if (!dir.exists(paste0(outputFolder, "/", collectionName, "_RES")))
      dir.create(paste0(outputFolder, "/", collectionName, "_RES"))
  }
  
  if( ("step_per_ranking" %in% names(config)) && (config[["step_per_ranking"]]))  {
    # add a column "rank_for_step" that holds 1-(rank/nrow(correlation_pairs_NPP))
    # it will be used as the "weight" to multiply the up step in the calculation of the ES
    nrow_correlation_pairs <- nrow(correlation_pairs_NPP)
    correlation_pairs_NPP[ ,rank_for_step:=1-.I/nrow_correlation_pairs]
    
    ranking_multiply <- config[["ranking_multiply"]]
    ES_values <- sapply(names(geneGroupsPairs), function(gene_set_name) { 
      RES_file_name <- ""
      if (write_RES_to_file)
        RES_file_name <- paste0(outputFolder, "/", collectionName, "_RES/RES_values_",
                                nppID, "NPP_", gene_set_name, config[["suffix"]], ".csv")
      GSEA.EnrichmentScore.per_ranking(gene_set_name, geneGroupsPairs[[gene_set_name]], 
                           correlation_pairs_NPP, nppID, write_RES_to_file,
                           RES_file_name, ranking_multiply)
    }, USE.NAMES = T)
  } else {
  ES_values <- sapply(names(geneGroupsPairs), function(gene_set_name) { 
    RES_file_name <- ""
    if (write_RES_to_file)
      RES_file_name <- paste0(outputFolder, "/", collectionName, "_RES/RES_values_", 
                              nppID, "NPP_", gene_set_name, config[["suffix"]], ".csv")
    GSEA.EnrichmentScore(gene_set_name, geneGroupsPairs[[gene_set_name]], 
                         correlation_pairs_NPP, nppID, write_RES_to_file, RES_file_name)
  }, USE.NAMES = T)
  }
  ES_dt <- data.table(gene_set=names(ES_values), 
                      ES=ES_values, 
                      gene_set_size=lengths(geneGroups[names(ES_values)]))
  fwrite(ES_dt, output_ES_file)
}
