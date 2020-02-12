
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


# for multiple combinations NPPs:
if "NPPs_from_file" in config and config["NPPs_from_file"]:
    with open(os.path.join(config["dataFolder"], config["combinations_list_file"])) as f:
        NPPs = f.read().splitlines()
else:
    NPPs = config["NPPs"]

include: "snakemake_ES_with_random.smk"

rule all_zscores:
    input:
        outputFolder + "/" + config["results_file_prefix"] + "_results.pdf"


rule mirror_random_ES_values:
    input:
        random_ES_file = outputFolder + "/" + randomFolder + "/ES_values_{NPP}NPP__{CORRELATION}_for_random_" +
                         collectionName + "_size_{geneSetSize}.csv"
    output:
        mirrored_random_ES_file = outputFolder + "/" + randomFolder + "/mirrored_ES_values_{NPP}NPP__{CORRELATION}_for_random_" +
                                  collectionName + "_size_{geneSetSize}.csv",
        mean_sd_file = temp(outputFolder + "/" + randomFolder + "/mirrored_mean_sd_{NPP}NPP__{CORRELATION}_for_random_" +
                            collectionName + "_size_{geneSetSize}.csv")
    script:
        "mirror_random_ES_snakemake.R"

rule merge_mean_sd_files:
    input:
        mean_sd_files = lambda wildcards: expand(outputFolder + "/" + randomFolder +
                                                 "/mirrored_mean_sd_{NPP}NPP__{CORRELATION}_for_random_" + collectionName + "_size_{geneSetSize}.csv",
                                                 NPP=wildcards.NPP, CORRELATION=wildcards.CORRELATION,
                                                 geneSetSize=config["collectionSizes"])
    output:
        filename = outputFolder + "/" + randomFolder + "/mirrored_mean_sd_{NPP}NPP__{CORRELATION}_for_random_" + collectionName + ".csv"
    run:
        with open(output.filename, "w+") as outfile:
            outfile.write("gene_set_size,random_ES_mean,random_ES_sd\n")
            for fname in input:
                with open(fname, 'r+') as infile:
                    shutil.copyfileobj(infile, outfile)


rule calculate_ES_zscores_for_random:
    input:
        NPP_ES_values_file = outputFolder + "/ES_values_{NPP}NPP_{CORRELATION}_" + collectionName + suffix + ".csv",
        mean_sd_file = rules.merge_mean_sd_files.output.filename
    output:
        outputFolder + "/ES_ZSCORES_{NPP}NPP_{CORRELATION}_for_" + randomFolder + ".csv"

    script:
        "calculate_ES_zscores_snakemake.R"

rule plot_ES_zscore_results:
    input:
        expand(outputFolder + "/ES_ZSCORES_{NPP}NPP_{CORRELATION}_for_" + randomFolder + ".csv",
               NPP=NPPs, CORRELATION=config["CORRELATIONs"])
    output:
        outputFolder + "/" + config["results_file_prefix"] + "_results.pdf"
    script:
        "plot_ES_zscores_results_snakemake.R"
