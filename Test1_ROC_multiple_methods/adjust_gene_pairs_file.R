library(data.table)

# USAGE:
# Rscript adjust_gene_pairs_file.R <linkages unsorted file> <sorted linkages file> <status letter: N or P>

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 3 arguments must be supplied (input gene pairs file, output file, Status)", call.=FALSE)
} 
 
gene_pairs_file <- args[1]
output_gene_pairs_file <- args[2]
status_val <- args[3]

gene_pairs <- fread(gene_pairs_file, header = F)

setnames(gene_pairs, "V1", "GeneA")
setnames(gene_pairs, "V2", "GeneB")

# add LinkID
gene_pairs[, LinkID := .I]
setcolorder(gene_pairs, c("LinkID", "GeneA", "GeneB"))

#  order GeneA, GeneB so that GeneA < GeneB alphabetically
unsorted_lines <- gene_pairs[GeneA > GeneB]
gene_pairs[LinkID %in% unsorted_lines$LinkID, c("GeneA", "GeneB") := list(GeneB, GeneA)]

# add column gene_pairs with GeneA:GeneB
gene_pairs[ ,gene_pair := paste0(GeneA, ":", GeneB)]

# add Status column (N for negative pairs, P for positive)
gene_pairs[ ,Status := status_val]

fwrite(gene_pairs, output_gene_pairs_file)




