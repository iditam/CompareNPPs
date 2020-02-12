#!/usr/bin/r
library(data.table)
library(yaml)
library(stringr)
library(plyr)

config <- yaml.load_file("config.yaml")

# transformation functions
one_divide_sqrt <- function(x) {
  1/sqrt(x)
}

one_divide_x <- function(x) {
  1/x
}

no_transform <- function(x) {
  x
}

create_combination_NPP <- function(threshold, transformation, divide_by_query, 
                                   colzscore, npp_filename, rowzscore="rowzscore") {
  
  if (str_detect(threshold, "late_flooring") || str_detect(threshold, "after_zscore")
      || str_detect(threshold, "below_") ) { # flooring will be made later
    matrix_file <- paste0(config$bitscoresFolder, "/", config$queryGenome, "_", config$NPP_resource,
                          "_BitScoreMatrix_no_threshold")
  } else {
    matrix_file <- paste0(config$bitscoresFolder, "/", config$queryGenome, "_", config$NPP_resource,
                          "_BitScoreMatrix_", threshold) # e.g. Homo_sapiens_UniProt_BitScoreMatrix_threshold_20.4
  }
  
  
  scores_matrix <- fread(matrix_file, header = T)
  matrix_rownames <- scores_matrix[ ,V1]
  # continue without rownames, just with scores
  scores_matrix <- scores_matrix[ ,2:ncol(scores_matrix)]
  
  # flooring to specific value, e.g.  "below_50_becomes_1"
  if (str_detect(threshold, "below_")) {
    vals <- str_match(threshold, "below_(.+)_becomes_(.+)")[2:3]
    threshold_val <- as.numeric(vals[1])
    new_val <- as.numeric(vals[2])
    scores_matrix[which(scores_matrix < threshold_val, arr.ind=T)] <- new_val
  }
  
  Protein2Genes <- fread(config$protein2gene_dict, header = F)
  setnames(Protein2Genes, c("V1", "V2"), c("Protein", "Gene"))
  
  # 1. convert from protein IDs to gene symbol
  gene_names <- mapvalues(matrix_rownames, Protein2Genes[ ,Protein], 
                          Protein2Genes[ ,Gene])
  
  # Preparation for threshold after z-score, if needed:
  # put NA in all values in P < threshold_val
  if (str_detect(threshold, "_after_zscore")) {
    threshold_val <- as.numeric(str_remove(str_remove(threshold, "_after_zscore"), "threshold_"))
    scores_matrix[scores_matrix < threshold_val] <- NA
  }
  
  # 2. divide values by query organism values - if needed, upon the value of "divide_by_query"
  if (divide_by_query == "divide") {
    queryOrgValuesReplicatedMatrix <- as.matrix(scores_matrix[ , get(config$queryGenomeId)]) %*% rep(1.0, ncol(scores_matrix))
    divided_scores_matrix <- scores_matrix/queryOrgValuesReplicatedMatrix
    # change all values > 1 to be 1
    divided_scores_matrix[divided_scores_matrix > 1] <- 1
  } else {
    divided_scores_matrix <- scores_matrix
  }
  
  # remove the query organism, not needed anymore 
  divided_scores_matrix <- divided_scores_matrix[ , -1]
  
  
  ## late flooring 
  if (str_detect(threshold, "late_flooring")) {
    threshold_val <- as.numeric(str_remove(threshold, "late_flooring_"))
    divided_scores_matrix[divided_scores_matrix < threshold_val] <- threshold_val
  }
  
  # 3. Transformation
  transformed_scores_matrix <- do.call(transformation, list(divided_scores_matrix))
  
  # 4. scale by columns - if needed, upon the value of "colzscore"
  if (colzscore == "colzscore") {
    col_scaled_matrix <- scale(transformed_scores_matrix)
  } else {
    col_scaled_matrix <- as.matrix(transformed_scores_matrix)
  }
  
  # if threshold after zscore needed - use the NA values we marked before,
  # and replace each NA value with the minimum value of the containing column
  if (str_detect(threshold, "_after_zscore")) {
    colmins <- integer(ncol(col_scaled_matrix))
    # num_of_NAs <- integer(ncol(col_scaled_matrix))
    for(i in 1:ncol(col_scaled_matrix)){
      colmins[i] <- min(col_scaled_matrix[,i], na.rm = TRUE)
      # num_of_NAs[i] <- length(col_scaled_matrix[is.na(col_scaled_matrix[,i]), i])
      col_scaled_matrix[is.na(col_scaled_matrix[,i]), i] <- min(col_scaled_matrix[,i], na.rm = TRUE)
    }
    # sum(num_of_NAs)
    # [1] 662500
  }
  
  # 5. scale by rows
  row_scaled_NPP <- t(scale(t(col_scaled_matrix)))
  # row_scaled_NPP <- col_scaled_matrix
  
  # add the gene column at the start of the table
  rownames(row_scaled_NPP) <- gene_names
  
  # fwrite() NPP
  if (!dir.exists(dirname(npp_filename))) {
    dir.create(dirname(npp_filename), showWarnings = FALSE)
  }
  
  npp_filename_taxid <- paste0(npp_filename, "_taxid")
  # npp_filename_taxid <- paste0(output_folder, "/", config$nppFilePrefix, "_", basename(output_folder), "_taxid_norowzscore")
  write.table(row_scaled_NPP, file=npp_filename_taxid, quote=F, col.names=NA, sep="\t")
  
  # convert to taxnames
  species_dict_file <- paste0(config$dataFolder, "/", config$species_list_file)
  species_dict <- fread(species_dict_file)
  colnames(row_scaled_NPP) <- mapvalues(colnames(row_scaled_NPP), species_dict[ ,taxid], species_dict[ ,taxname], warn_missing = F)
  # npp_filename <- paste0(output_folder, "/", config$nppFilePrefix, "_", basename(output_folder), "_norowzscore")
  
  write.table(row_scaled_NPP, file=npp_filename, quote=F, col.names=NA, sep="\t")
  
  
}
