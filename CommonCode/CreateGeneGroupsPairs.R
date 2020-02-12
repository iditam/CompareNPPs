#!/usr/bin/Rscript
#SBATCH --mem-per-cpu=2000m
#SBATCH -t1-0
#SBATCH -c16
#SBATCH --output=logs/CompareNPPs.%j.log
#SBATCH --error=logs/CompareNPPs.%j.log

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

sorted_gene_pairs <- function(geneGroup) {
  all_pairs <-  combn(geneGroup, 2, alphabetical_paste_vec)
  sort(all_pairs)
}

# Now create a vector with all gene pairs that are found in any gene set, in same format: "gene1:gene2",
# without repeats
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl, c("alphabetical_paste", "alphabetical_paste_vec", "sorted_gene_pairs")) # functions/variables
geneGroupsPairs <- parLapply(cl, geneGroups, sorted_gene_pairs)
stopCluster(cl)
saveRDS(geneGroupsPairs, paste0(geneGroupsFolder, "/geneGroupsPairs", "_", config[["collectionName"]], ".rds"))

# 
# # unique and sort gene groups pairs ----
# 
# # Run this part with a lot of memory
# geneGroupsPairs <- readRDS(paste(outputRdsFolder, "/geneGroupsPairs", "_", geneGroupsCollection, ".rds", sep = ""))
# allGeneGroupsPairs <- sort(unique(unlist(geneGroupsPairs)))
# saveRDS(allGeneGroupsPairs, paste(outputRdsFolder, "/allGeneGroupsPairs", "_", geneGroupsCollection, ".rds", sep = ""))
#  
