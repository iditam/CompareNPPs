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


include: "snakemake_correlations.smk"

rule all:
    input:
        ROC_filename = os.path.join(config["outputFolder"], config["ROC_file_prefix"] + "_" + config["output_files_id"] + ".pdf"),
        PR_filename = os.path.join(config["outputFolder"], config["PR_file_prefix"] + "_" + config["output_files_id"] + ".pdf"),
        AUCs_filename = os.path.join(config["outputFolder"], config["AUCs_file_prefix"] + "_" + config["output_files_id"] + ".csv")


rule create_linkages_with_npp_correlations:
    input:
        linkage_file = os.path.join(config["linkagesFolder"], "{LINKAGE_PREFIX}" + ".csv"),
        correlation_file = os.path.join(dataFolder, "{NPP}", "{CORRELATION}_correlation_pairs_{NPP}NPP.csv")
    output:
        filename = os.path.join(config["linkagesFolder"],
                                "{LINKAGE_PREFIX}" + "-" + "{NPP}" + "-" + "{CORRELATION}" + ".csv")
    script:
        "create_npp_linkages_snakemake.R"

rule merge_positive_linkages:
    input:
        linkage_file = os.path.join(config["linkagesFolder"], config["positiveLinkagesPrefix"] + ".csv"),
        npp_linkage_files = expand(os.path.join(config["linkagesFolder"],
                                   config["positiveLinkagesPrefix"] + "-" + "{NPP}" + "-" + "{CORRELATION}" + ".csv"),
                                   NPP=NPPs, CORRELATION=config["CORRELATIONs"])
    output:
        filename = os.path.join(config["linkagesFolder"],
                                config["positiveLinkagesPrefix"] + "_" + config["output_files_id"] + ".csv")
    script:
        "merge_npps_linkages_snakemake.R"

rule merge_negative_linkages:
    input:
        linkage_file = os.path.join(config["linkagesFolder"], config["negativeLinkagesPrefix"] + ".csv"),
        npp_linkage_files = expand(os.path.join(config["linkagesFolder"],
                                   config["negativeLinkagesPrefix"] + "-" + "{NPP}" + "-" + "{CORRELATION}" + ".csv"),
                                   NPP=NPPs, CORRELATION=config["CORRELATIONs"])
    output:
        filename = os.path.join(config["linkagesFolder"],
                                config["negativeLinkagesPrefix"] + "_" + config["output_files_id"] + ".csv")
    script:
        "merge_npps_linkages_snakemake.R"

rule ROC_NPP_methods:
    input:
        positive_linkages = os.path.join(config["linkagesFolder"],
                                         config["positiveLinkagesPrefix"] + "_" + config["output_files_id"] + ".csv"),
        negative_linkages = os.path.join(config["linkagesFolder"],
                                         config["negativeLinkagesPrefix"] + "_" + config["output_files_id"] + ".csv")
    output:
        ROC_filename = os.path.join(config["outputFolder"], config["ROC_file_prefix"] + "_" + config["output_files_id"] + ".pdf"),
        PR_filename = os.path.join(config["outputFolder"], config["PR_file_prefix"] + "_" + config["output_files_id"] + ".pdf"),
        AUCs_filename = os.path.join(config["outputFolder"], config["AUCs_file_prefix"] + "_" + config["output_files_id"] + ".csv")
    script:
        "ROC_NPP_methods_snakemake.R"
