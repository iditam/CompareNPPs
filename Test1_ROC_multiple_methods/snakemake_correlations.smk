##################################################################################
# before running this snakemake run manualy "adjust_linkage_files_to_NPP_genes.R"
# to add gene symbols GeneA,GeneB to positive & negative linkages files,
# to sort them so that GeneA < GeneB, and to adjust gene symbols to NPP
# (convert protein entrez ids appropriately)
##################################################################################

import os
import shutil

outputFolder = config["outputFolder"]
dataFolder = config["dataFolder"]
nppFilePrefix = config["nppFilePrefix"]
numOfGenes = config["numOfGenes"]
if "numOfBins" in config:
    numOfBins = config["numOfBins"]
else:
    numOfBins = 500

# create list of pairs of indices for the rule create_partial_correlation_pairs
numOfCorrelations = numOfGenes*(numOfGenes-1)//2
bin_size = numOfCorrelations // numOfBins
maxIndex = (numOfCorrelations // bin_size) + (numOfCorrelations % bin_size > 0)

index_pairs = []
for i in range(1, maxIndex+1):
    # print "range ",i, ":"
    min_index = (i-1)*bin_size + 1
    max_index = i*bin_size if (i*bin_size < numOfCorrelations) else numOfCorrelations
    index_pairs.append((min_index, max_index))

# for multiple combinations NPPs
if "NPPs_from_file" in config and config["NPPs_from_file"]:
    with open(os.path.join(config["dataFolder"], config["combinations_list_file"])) as f:
        NPPs = f.read().splitlines()
else:
    NPPs = config["NPPs"]


rule all_correlations:
    input:
        expand(os.path.join(dataFolder, "{NPP}", "{CORRELATION}_correlation_pairs_{NPP}NPP.csv"),
               NPP=NPPs, CORRELATION=config["CORRELATIONs"])

rule create_partial_correlation_pairs:
    input:
        filename = os.path.join(dataFolder, "{NPP}", config["nppFilePrefix"] + "_" + "{NPP}")
    output:
        filename = temp(os.path.join(dataFolder, "{NPP}", "{CORRELATION}_correlation_pairs_{NPP}NPP_{partial_index}.csv"))
    params:
        current_range = lambda wildcards: index_pairs[int(wildcards.partial_index)]
    script:
        "create_partial_correlation_pairs.py"
        # with open(output.filename,"w+") as text_file:
        #     text_file.write("%s" % wildcards.index_pair)
        # Rscript create_partial_correlation_pairs.R - -npp {wildcards.NPP} - -cor {wildcards.CORRELATION} - -range {wildcards.index_pair}


rule merge_correlation_pairs_files:
    input:
        lambda wildcards: expand(os.path.join(dataFolder, "{NPP}", "{CORRELATION}_correlation_pairs_{NPP}NPP_{partial_index}.csv"),
                                 NPP=wildcards.NPP, CORRELATION=wildcards.CORRELATION,
                                 partial_index=range(len(index_pairs)))

    output:
        filename = os.path.join(dataFolder, "{NPP}", "{CORRELATION}_correlation_pairs_{NPP}NPP.csv")
    run:
        with open(output.filename, "w+") as outfile:
            outfile.write("correlation,gene_pair\n")
            for fname in input:
                with open(fname, 'r+') as infile:
                    shutil.copyfileobj(infile, outfile)
