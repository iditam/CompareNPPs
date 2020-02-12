# library(yaml)
library(data.table)
library(doParallel)

# config <- yaml.load_file("config_with_random.yaml")

num_of_cores <- detectCores()
registerDoParallel(cores=num_of_cores)

set.seed(123)

create_random_pairs <- function(gene_set_size, all_genes){
  random_genes <- sample(all_genes, gene_set_size, replace = F)
  dt <- combn(random_genes, 2)
  # order so that GeneA < GeneB alphabetically
  dt <- transpose(data.table(apply(dt, 2, sort)))
  dt[ ,gene_pair := paste0(V1,":",V2)]
  dt[,gene_pair]
}

create_random_group_pairs <- function(gene_set_size, all_genes_file, output_pairs_file, config) {
  numOfGroupsPerSize <- config$randomTimes
  
  all_genes <- fread(file = all_genes_file, header = F)
  all_genes <- all_genes[ ,V1] # character vector with gene names

  pairs_mat <- foreach(i = 1:numOfGroupsPerSize, .combine=rbind, 
                       .export = c("gene_set_size", "all_genes")) %dopar% {
                         random_group_name <- paste0("RANDOM_GROUP_", i)
                         c(random_group_name, create_random_pairs(gene_set_size, all_genes))
                       }
  fwrite(as.data.table(pairs_mat), file = output_pairs_file, col.names = F)
}



