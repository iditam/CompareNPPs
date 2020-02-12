#!/usr/bin/Rscript
#SBATCH --mem-per-cpu=2000m
#SBATCH -t02:59:00
#SBATCH -c16
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log

# create gene groups pairs ----
# Run this script with many CPUs

library(yaml)
library(parallel)

config <- yaml.load_file("config.yaml")

source("../CommonCode/ReadGmt.R")

# Take the collection details from collectionsFile, upon collectionID
geneGroupsFolder = file.path(config[["geneGroupsRootFolder"]], config[["geneCollectionResource"]])
geneGroups <- get_genes_groups_list(geneGroupsFolder,config[["collectionName"]])

alphabetical_paste <- function(a,b) {
  ifelse(a < b, paste(a, b, sep = ":"), paste(b, a, sep = ":"))
}
alphabetical_paste_vec <- function(vec) {
  a <- vec[1]
  b <- vec[2]
  alphabetical_paste(a,b)
}

sorted_gene_pairs_string <- function(geneGroup) {
  all_pairs <-  combn(geneGroup, 2, alphabetical_paste_vec)
  paste(sort(as.character(all_pairs)), collapse= ",")
}

# Now create a vector with all gene pairs that are found in any gene set, in same format: "gene1:gene2",
# without repeats
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl, c("alphabetical_paste", "alphabetical_paste_vec", "sorted_gene_pairs_string")) # functions/variables
geneGroupsPairsStrings <- parSapply(cl, geneGroups, sorted_gene_pairs_string)
stopCluster(cl)
pairs_file <- paste0(geneGroupsFolder, "/geneGroupsPairs", "_", config[["collectionName"]], ".gmt")
lines <- paste(names(geneGroupsPairsStrings), geneGroupsPairsStrings, sep = ",")
writeLines(lines, con = pairs_file)
# lapply(names(geneGroupsPairs), function(group_name) { write(c(group_name, geneGroupsPairs[[group_name]]),
#                                                             file = pairs_file, sep = ",", append = T)})


# 
# # unique and sort gene groups pairs ----
# 
# # Run this part with a lot of memory
# geneGroupsPairs <- readRDS(paste(outputRdsFolder, "/geneGroupsPairs", "_", geneGroupsCollection, ".rds", sep = ""))
# allGeneGroupsPairs <- sort(unique(unlist(geneGroupsPairs)))
# saveRDS(allGeneGroupsPairs, paste(outputRdsFolder, "/allGeneGroupsPairs", "_", geneGroupsCollection, ".rds", sep = ""))
#  
