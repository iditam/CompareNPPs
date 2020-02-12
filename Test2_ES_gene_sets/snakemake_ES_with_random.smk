
import shutil


outputFolder = config["outputFolder"]
dataFolder = config["dataFolder"]
allGenesFile = config["allGenesFile"]
geneGroupsRootFolder = config["geneGroupsRootFolder"]
geneCollectionResource = config["geneCollectionResource"]
collectionName = config["collectionName"]
suffix = config["suffix"]
write_RES_to_file = config["write_RES_to_file"]
randomTimes = config["randomTimes"]
randomBinSize = config["randomBinSize"]

randomFolder = collectionName + "_random_" + str(randomTimes) + "_ES" + suffix

# randomTimes = 1000
# randomBinSize = 10

maxIndex = (randomTimes // randomBinSize) + (randomTimes % randomBinSize > 0)

index_pairs = []
for i in range(1, maxIndex+1):
    # print "range ",i, ":"
    min_index = (i-1)*randomBinSize + 1
    max_index = i*randomBinSize if (i*randomBinSize < randomTimes) else randomTimes
    index_pairs.append((min_index, max_index))

# for multiple combinations NPPs:
if "NPPs_from_file" in config and config["NPPs_from_file"]:
    with open(os.path.join(config["dataFolder"], config["combinations_list_file"])) as f:
        NPPs = f.read().splitlines()
else:
    NPPs = config["NPPs"]


rule all:
    input:
        expand(outputFolder + "/ES_{NPP}NPP_{CORRELATION}_with_" + randomFolder + ".pdf",
               NPP=NPPs, CORRELATION=config["CORRELATIONs"])


rule calculate_ES:
    input:
        correlation_file = dataFolder + "/{NPP}/{CORRELATION}_correlation_pairs_{NPP}NPP.csv",
        gene_group_pairs = geneGroupsRootFolder + "/" + geneCollectionResource + "/geneGroupsPairs_" + collectionName + ".rds"
    output:
        outputFolder + "/ES_values_{NPP}NPP_{CORRELATION}_" + collectionName + suffix + ".csv"
    script:
        "calculate_ES_snakemake.R"


rule order_and_manipulate_correlation_pairs:
    input:
        correlation_file = dataFolder + "/{NPP}/{CORRELATION}_correlation_pairs_{NPP}NPP.csv"
    output:
        dataFolder + "/{NPP}/{CORRELATION}_correlation_pairs_{NPP}NPP_ordered" + suffix + ".csv"
    script:
        "order_and_manipulate_correlation_pairs_snakemake.R"


rule create_random_group_pairs:
    input:
        all_genes_file = dataFolder + "/" + allGenesFile
    output:
        filename = protected(geneGroupsRootFolder + "/" + geneCollectionResource + "/" + collectionName + "_random_pairs_" + str(randomTimes) + "/random_pairs_size_{geneSetSize}.csv")
    script:
        "create_random_group_pairs_snakemake.R"


rule partial_calculate_ES_for_random_groups:
    input:
        correlation_file = rules.order_and_manipulate_correlation_pairs.output,
        random_pairs_file = rules.create_random_group_pairs.output.filename
    output:
        temp(outputFolder + "/" + randomFolder + "/partial_ES_values_{NPP}NPP__{CORRELATION}_for_random_" + collectionName + "_size_{geneSetSize}_{partial_index}.csv")
    # params:
    #     current_index = lambda wildcards: int(wildcards.partial_index)
    params:
        current_range = lambda wildcards: index_pairs[int(wildcards.partial_index)]
    script:
        "calculate_ES_for_random_groups_partial_snakemake.R"

rule merge_ES_values_files:
    input:
        lambda wildcards: expand(outputFolder + "/" + randomFolder +
                                 "/partial_ES_values_{NPP}NPP__{CORRELATION}_for_random_" + collectionName +
                                 "_size_{geneSetSize}_{partial_index}.csv",
                                 NPP=wildcards.NPP, CORRELATION=wildcards.CORRELATION,
                                 geneSetSize=wildcards.geneSetSize,
                                 partial_index=range(len(index_pairs)))
                                 # partial_index=range(1, randomTimes+1))

    output:
        filename = outputFolder + "/" + randomFolder +
            "/ES_values_{NPP}NPP__{CORRELATION}_for_random_" + collectionName + "_size_{geneSetSize}.csv"
    run:
        with open(output.filename, "w+") as outfile:
            for fname in input:
                with open(fname, 'r+') as infile:
                    shutil.copyfileobj(infile, outfile)

rule plot_ES_results_with_random:
    input:
        NPP_ES_values_file = rules.calculate_ES.output,
        random_ES_values_files = lambda wildcards: expand(outputFolder + "/" + randomFolder + "/ES_values_{NPP}NPP__{CORRELATION}_for_random_" + collectionName + "_size_{geneSetSize}.csv",
                                                          NPP=wildcards.NPP, CORRELATION=wildcards.CORRELATION,
                                                          geneSetSize=config["collectionSizes"])
    output:
        outputFolder + "/ES_{NPP}NPP_{CORRELATION}_with_" + randomFolder + ".pdf"

    script:
        "plot_ES_with_random_snakemake.R"
