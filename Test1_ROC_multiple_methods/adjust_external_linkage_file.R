##########################################################
# This script processes the files of positive linkages and negative linkages (with given names)
# These files are csv with just pairs - each line is of the form: Gene1,Gene2
# and adjusts them for usage in the snakemake pipeline:
# 1. order each line GeneA,GeneB so that GeneA < GeneB alphabetically
# 2. add a column "gene_pair" with <GeneA>:<GeneB> (just for convenience for further searches)
# 3. remove all duplicate linkages (rows with same "gene_pair")
##### TBD: ####
# save the results files upon the values <positiveLinkagesFile>,<negativeLinkages123File> in config.yaml 
##########################################################

library(yaml)
library(data.table)
library(plyr)
library(stringr)

config <- yaml.load_file("config.yaml")

nppID = "UniProt_062018"
protein2gene_file = "Homo_sapiens_uniprot_id2gene_062018.csv"
original_positive_linkages_file = "positive_linkages_original.csv"
original_negative123_linkages_file = "negative_linkages_123_original.csv"

pos_linkages <- fread(file = file.path(config$linkagesFolder, original_positive_linkages_file))
neg_linkages_123 <- fread(file = file.path(config$linkagesFolder, original_negative123_linkages_file))


# 1. for the rows that have gene symbols that do not appear in NPP,
# fetches a mapping of the entrez id's

# check how many gene pairs do not exist in NPP
UniProtNPP_protein2genes <- fread(file = file.path(config$dataFolder, nppID, protein2gene_file), header = F)
setnames(UniProtNPP_protein2genes, "V1", "uniprot_ac")
setnames(UniProtNPP_protein2genes, "V2", "gene")
UniProtNPP_genes = UniProtNPP_protein2genes$gene
length(UniProtNPP_genes)

length(which( (!pos_linkages$GeneA %in% UniProtNPP_genes) 
              | (!pos_linkages$GeneB %in% UniProtNPP_genes)))
# 2058 out of 57114 pairs
# 2058/57114
# [1] 0.0360332

length(which( (!neg_linkages_123$GeneA %in% UniProtNPP_genes)
              | (!neg_linkages_123$GeneB %in% UniProtNPP_genes)))
# 37856 out of 571140 pairs
# 37856/571140
# [1] 0.06628147

unconverted_entrez_ids <- unique(c(pos_linkages[!GeneA %in% UniProtNPP_genes, ProteinA], 
                                   pos_linkages[!GeneB %in% UniProtNPP_genes, ProteinB],
                                   neg_linkages_123[!GeneA %in% UniProtNPP_genes, ProteinA], 
                                   neg_linkages_123[!GeneB %in% UniProtNPP_genes, ProteinB]))
unconverted_entrez_ids <- str_remove(unconverted_entrez_ids,"hsa:")
writeLines(unconverted_entrez_ids, file.path(config$linkagesFolder, "unconverted_entrez_ids.csv"))

### manually fetched mapping by uploading unconverted_entrez_ids.csv to https://www.uniprot.org/uploadlists/
### The mapping was saved to the file unconverted_entrez_ids_mapping.csv

# build a dictionary protein_id <--> gene_in_npp based on unconverted_entrez_ids_mapping.csv
missing_proteins_mapping <- data.table(entrez_id = unconverted_entrez_ids, 
                                       uniprot_ac = character(length(unconverted_entrez_ids)),
                                       gene_in_npp = character(length(unconverted_entrez_ids)))

mapping_from_uniprot <- fread(file = file.path(config$linkagesFolder, "unconverted_entrez_ids_mapping.csv"))

setnames(mapping_from_uniprot, colnames(mapping_from_uniprot)[1], "entrez_id")
setnames(mapping_from_uniprot, "Entry", "uniprot_ac")

found_uniprot_ac_mapping <- mapping_from_uniprot[uniprot_ac %in% UniProtNPP_protein2genes$uniprot_ac, .(entrez_id, uniprot_ac)]

missing_proteins_mapping$uniprot_ac = mapvalues(missing_proteins_mapping$entrez_id,
                                                 found_uniprot_ac_mapping$entrez_id,
                                                 found_uniprot_ac_mapping$uniprot_ac, warn_missing = F)

missing_proteins_mapping$gene_in_npp = mapvalues(missing_proteins_mapping$uniprot_ac,
                                                 UniProtNPP_protein2genes$uniprot_ac,
                                                 UniProtNPP_protein2genes$gene, warn_missing = F)

# Now deal manually with the entrez_id's that could not find a mapping to uniprot_ac that appear in UniProtNPP_protein2genes
# They have identical values in entrez_id,uniprot_ac,gene_in_npp as the 'mapvalues' did not find a convertion
missing_proteins_mapping[entrez_id==uniprot_ac]
# entrez_id uniprot_ac gene_in_npp
# 1:      6171       6171        6171

# This entrez_id has the genes: RPL41 hCG_2039605
 "RPL41" %in% UniProtNPP_genes
 # [1] FALSE
 "hCG_2039605" %in% UniProtNPP_genes
 # [1] FALSE
 
 # not found, it's just 1 protein, we will remove it from missing_proteins_mapping 
 missing_proteins_mapping <- missing_proteins_mapping[entrez_id != 6171]
 
 # Now convert the gene symbols of missing_proteins_mapping in pos_linkages, neg_linkages_123 
 # pos_linkages_before_convert <- fread(file.path(config$linkagesFolder, config$positiveLinkagesFile))
 pos_linkages_before_convert <- copy(pos_linkages)
 
 get_npp_gene <- function(entrezid, original_gene) {
   gene_in_npp <- missing_proteins_mapping[entrez_id==entrezid, gene_in_npp][1]
   if(is.na(gene_in_npp))
     gene_in_npp <- original_gene
   return(gene_in_npp)
 }

pos_linkages[!GeneA %in% UniProtNPP_genes, GeneA := mapply(get_npp_gene, str_remove(ProteinA, "hsa:"), GeneA)]
unique(pos_linkages[!GeneA %in% UniProtNPP_genes, ProteinA])
# [1] "hsa:6171" - indeed this is the protein we did not manage to map..
pos_linkages[!GeneB %in% UniProtNPP_genes, GeneB := mapply(get_npp_gene, str_remove(ProteinB, "hsa:"), GeneB)]
unique(pos_linkages[!GeneB %in% UniProtNPP_genes, ProteinB])
# [1] "hsa:6171" - indeed this is the protein we did not manage to map..
length(which( (!pos_linkages$GeneA %in% UniProtNPP_genes) 
              | (!pos_linkages$GeneB %in% UniProtNPP_genes)))
# 126 out of 57114 pairs
# 126/57114
# [1] 0.002206114

neg_linkages_123[!GeneA %in% UniProtNPP_genes, GeneA := mapply(get_npp_gene, str_remove(ProteinA, "hsa:"), GeneA)]
unique(neg_linkages_123[!GeneA %in% UniProtNPP_genes, ProteinA])
# character(0)
neg_linkages_123[!GeneB %in% UniProtNPP_genes, GeneB := mapply(get_npp_gene, str_remove(ProteinB, "hsa:"), GeneB)]
unique(neg_linkages_123[!GeneB %in% UniProtNPP_genes, ProteinB])
# character(0)
length(which( (!neg_linkages_123$GeneA %in% UniProtNPP_genes) 
              | (!neg_linkages_123$GeneB %in% UniProtNPP_genes)))
# [1] 0


# 2. order ProteinA, ProteinB, GeneA, GeneB so that GeneA < GeneB alphabetically
unsorted_lines <- pos_linkages[GeneA > GeneB]
pos_linkages[LinkID %in% unsorted_lines$LinkID, c("GeneA", "GeneB", "ProteinA", "ProteinB") := 
               list(GeneB, GeneA, ProteinB, ProteinA)]
# sanity check:
a <- pos_linkages[LinkID %in% unsorted_lines$LinkID, c(1,3,2,16,15)]
b <- unsorted_lines[ ,c(1,2,3,15,16)]
colnames(a) <- as.character(1:5)
colnames(b) <- as.character(1:5)
all.equal(a,b)
# [1] TRUE - OK


# negative linkages

unsorted_lines <- neg_linkages_123[GeneA > GeneB]
neg_linkages_123[LinkID %in% unsorted_lines$LinkID, c("GeneA", "GeneB", "ProteinA", "ProteinB") := 
                   list(GeneB, GeneA, ProteinB, ProteinA)]
# sanity check:
a <- neg_linkages_123[LinkID %in% unsorted_lines$LinkID, c(1,3,2,16,15)]
b <- unsorted_lines[ ,c(1,2,3,15,16)]
colnames(a) <- as.character(1:5)
colnames(b) <- as.character(1:5)
all.equal(a,b)
# [1] TRUE

# 3. add a column "gene_pair" with <GeneA>:<GeneB>
pos_linkages[ ,gene_pair := paste0(GeneA, ":", GeneB)]
neg_linkages_123[ ,gene_pair := paste0(GeneA, ":", GeneB)]

# keep a copy of linkages before unique
fwrite(pos_linkages, file =  file.path(config$linkagesFolder, paste0(config$positiveLinkagesFile, "_non_unique")))
fwrite(neg_linkages_123, file =  file.path(config$linkagesFolder, paste0(config$negativeLinkages123File, "_non_unique")))

# 4. remove all duplicate linkages (rows with same "gene_pair")
pos_linkages <- unique(pos_linkages,by = "gene_pair")
fwrite(pos_linkages, file =  file.path(config$linkagesFolder, config$positiveLinkagesFile))
neg_linkages_123 <- unique(neg_linkages_123,by = "gene_pair")
fwrite(neg_linkages_123, file =  file.path(config$linkagesFolder, config$negativeLinkages123File))



