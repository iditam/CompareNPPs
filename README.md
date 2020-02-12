How to use the code for testing a version of NPP

1. in folder NPP_data: create a sub-folder and put the npp table file inside it

2. In folder Gene_sets: create a sub-folder that contains a "Resource" of gene sets (e.g KEGG, MSigBB..).

3. In folder Gene_sets\your_resource put:
  a. a .rds file that contains a list of named character vectors, each representing a gene set, with it's name.
  b. a .csv file that lists all the gene set collections (see KeggFromAPI\kegg_files.csv for example)

4. copy all files from the CommonCode dataFolder

5. To create the geneGroupsPairs file (e.g. geneGroupsPairs_kegg_pathway_genes.rds) :
edit CommonCode\config.yaml
run:
> cd CommonCode
> Rscript CreateGeneGroupsPairs.R

Test1 (ROC)
----------
1. Copy all files from Test1_ROC_multiple_methods and create sub-folders "logs", "Results"

2. install snakemake:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

3. cd Test1_ROC_multiple_methods
edit your config.yaml
edit snakemake_icore.json
run the test:
bash run_snakemake_ROC_multiple.sh
to track the running, check the files:
nohup.out
logs\*.log

The results will appear in Test1_ROC_multiple_methods\Results

Test2 (ES)
----------
1. Copy all files from Test2_ES_gene_sets and create sub-folders "logs", "Results"

2. install snakemake:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

3. cd Test2_ES_gene_sets
edit your config.yaml
edit snakemake_ES_with_random_icore.json
run the test:
bash run_snakemake_zscores.sh
to track the running, check the files:
zscores_snakemake.out
logs\*.log

The results will appear in Test2_ES_gene_sets\Results

