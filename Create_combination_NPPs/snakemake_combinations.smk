

rule all:
    input:
        config["dataFolder"] + "/" + config["combinations_list_file"],
        config["dataFolder"] + "/" + config["combinations_details_file"]

rule create_combination_NPP:
    input:
    output:
        npp_filename = config["dataFolder"] + "/" + config["NPP_resource"] +
                               "_{threshold}__{transformation}__{divide_by_query}__{colzscore}/" +
                               config["queryGenome"] + "_NPP_" + config["NPP_resource"] +
                                                      "_{threshold}__{transformation}__{divide_by_query}__{colzscore}"
    script:
        "create_combination_NPP_snakemake.R"

rule create_combinations_list:
    input:
        expand(config["dataFolder"] + "/" + config["NPP_resource"] +
               "_{threshold}__{transformation}__{divide_by_query}__{colzscore}/" +
                config["queryGenome"] + "_NPP_" + config["NPP_resource"] +
                                       "_{threshold}__{transformation}__{divide_by_query}__{colzscore}",
               threshold=config["thresholds"], transformation=config["transformations"],
               divide_by_query=config["divide_by_query_list"],
               colzscore=config["colzscore_list"])
    output:
        filename = temp(config["dataFolder"] + "/" + config["temp_combinations_list_file"])
    run:
        with open(output.filename, "w+") as outfile:
            for npp_path in input:
                # npp_file = "/".join(npp_path.strip("/").split('/')[1:])  # remove the beginning: config["dataFolder"] + "/"
                npp_file = os.path.basename(npp_path)
                outfile.write(npp_file + "\n")

rule update_combinations_details:
    input:
        filename = config["dataFolder"] + "/" + config["temp_combinations_list_file"]
    output:
        list_file = config["dataFolder"] + "/" + config["combinations_list_file"],
        details_file = config["dataFolder"] + "/" + config["combinations_details_file"]
    script:
        "update_combination_details_snakemake.R"
