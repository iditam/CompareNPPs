

# same args as lapply.  FUN must return named components. 
lapplyWithRename <- function(...) { 
  x <- lapply(...) 
  names(x) <- sapply(x, names) 
  lapply(x, function(x) { names(x) <- NULL; x }) 
} 

set_line_names <- function(x) {
  splittedLine <- strsplit(x, "\t")[[1]]
  structure(x, .Names = splittedLine[1]) 
}

split_gene_groups <- function(x) {
  splittedLine <- strsplit(x, "\t")[[1]]
  x <- splittedLine[c(3:length(splittedLine))]
  x
}

# Reads a 'gmt' format file, each line representing a gene group, with fields separated with '\t' 
# 1st field: group name. from 3st field until end of line: genes
# Creates a list of character vectors, each vector is a gene group, 
# the names in the list are the group names
# Saves a 'rds' file of the list, and return it as return value

get_genes_groups_list<-function(gmtFolder, collectionName){
  rdsFile <- paste0(gmtFolder, "/", collectionName, ".rds")
  if (file.exists(rdsFile))
    geneGroupsList <- readRDS(rdsFile)
  else {
    genesGmtFile <- paste(gmtFolder, "/", collectionName, ".gmt", sep = "")
    gmtFileLines<-scan(genesGmtFile, what=character(), sep="\n", multi.line=T, quiet = T)
    gmtLines <- lapplyWithRename(gmtFileLines, set_line_names) 
    geneGroupsList <- sapply(gmtLines, split_gene_groups, USE.NAMES = T)
    saveRDS(geneGroupsList, rdsFile)
  }
  geneGroupsList
}

get_all_lists <- function(gmtFolder, csvFile) {
  collectionDetails <- read.csv(paste(gmtFolder, "/", csvFile, sep = ""), as.is = T)
  all_lists <- sapply(collectionDetails$Name, function(n){get_genes_groups_list(gmtFolder, n)})
}

# gmtFolder <- "/Users/yuval/Google Drive/Data/MSigDB"
# collectionName <- "c1.all.v6.0.symbols"
# l <- get_genes_groups_list(gmtFolder, collectionName)
# l[["chr13q11"]]
