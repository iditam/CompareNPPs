
- create a folder for the new resource, e.g: ../Gene_sets/KeggFromAPI/
- if you have  .rds file[s] for the gene sets colletions, put them in that folder
  e.g. 
  kegg_pathways <- readRDS("../Gene_sets/KeggFromAPI/kegg_pathway_genes.rds")
  
- Create in that folder a .csv file with same format as ../Gene_sets/MSigDB/gmt_files.csv
- for each collection in the new resource, create a .rds of all gene pairs, to be read like this:
geneGroupsPairs <- readRDS(paste0(geneGroupsFolder, "/geneGroupsPairs_", collectionName, ".rds"))
Use the script CreateGeneGroupsPairs.sh (uses CreateGeneGroupsPairs.R)

