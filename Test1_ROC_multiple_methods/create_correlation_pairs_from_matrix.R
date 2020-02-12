library(data.table)


filename <- "BPP_hamming_Eval.csv"
similarity_name <- "hamming"
NPP_id <- "BPP_Eval"
dataFolder <- "../NPP_data"

cor_list_file <- paste0(dataFolder, "/", NPP_id, "/", similarity_name, 
                          "_correlation_pairs_", NPP_id, "NPP_datatable.csv")

mat_dt <- fread(paste0(dataFolder, "/", filename), header = T)
gene_names <- colnames(mat_dt)[-1]

genes_num <- length(gene_names)

# file_counter <- 1
# bin_size <- 4 # 10000  
# 
# correlations <- numeric()
# gene_pairs <- character()
# sapply(1:(genes_num-1), function(i) {
#    sapply((i+1):genes_num, function(j) {
# # sapply(1:5, function(i) {
# #   sapply((i+1):6, function(j) {
#     correlations <<- c(correlations, unlist(mat_dt[i, j+1, with=F]))
#     if (gene_names[i] < gene_names[j]) {
#       gene_pairs <<- c(gene_pairs, paste0(gene_names[i], ":", gene_names[j]))
#     } else {
#       gene_pairs <<- c(gene_pairs, paste0(gene_names[j], ":", gene_names[i]))
#     }
#   })
#   # flush to file
#   if (length(correlations) > bin_size) {
#     dt <- data.table(correlation=correlations,
#                      gene_pair=gene_pairs)
#     fwrite(dt, file = paste0(cor_list_prefix, file_counter, ".csv"))
#     correlations <<- numeric()
#     gene_pairs <<- character()
#     file_counter <<- file_counter + 1
#   }
# })

mat <- as.matrix(mat_dt[ ,2:ncol(mat_dt)])

mat[lower.tri(mat,diag=TRUE)]=NA
mat_table <- as.table(mat)
rm(mat)

rm(mat_dt)

rownames(mat_table) <- colnames(mat_table)
dt <- as.data.table(mat_table)
dt=na.omit(dt)

alphabetic_pair <- function(str1, str2) {
  if (str1 < str2) {
    return(paste0(str1, ":", str2))
  } else {
    return(paste0(str2, ":", str1))
  }
}

fwrite(dt, cor_list_file)

# fix file with:
# awk -F ',' '{if (NR==1) print "correlation,gene_pair"; else if ($1 < $2) print $3 "," $1 ":" $2; else print $3 "," $2 ":" $1; }' dist_correlation_pairs_PPPNPP_datatable.csv > dist_correlation_pairs_PPPNPP.csv
# awk -F ',' '{if (NR==1) print "correlation,gene_pair"; else if ($1 < $2) print $3 "," $1 ":" $2; else print $3 "," $2 ":" $1; }' hamming_correlation_pairs_BPP_EvalNPP_datatable.csv > hamming_correlation_pairs_BPP_EvalNPP.csv

# save.image("create_list.RData")
