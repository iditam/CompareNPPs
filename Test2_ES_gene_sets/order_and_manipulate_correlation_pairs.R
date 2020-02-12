# library(yaml)
library(data.table)
library(stringr)

# config <- yaml.load_file("config_with_random.yaml")

order_and_manipulate_correlation_pairs <- function(correlation_file, output_correlation_file, config) {
  # # TBD: remove
  # correlation_file <- "../NPP_data/UniProt_log2_before_divide/pearson_correlation_pairs_UniProt_log2_before_divideNPP.csv"

  # load correlation pairs
  correlation_pairs_NPP <- fread(input = correlation_file, verbose = F, showProgress = F)

  # order correlation pairs by correlation (descending)
  setorder(correlation_pairs_NPP, -correlation)
  
  if (config$suffix == "_zeroed") {
    # zero negative correlations
    correlation_pairs_NPP[correlation < 0, correlation:=0]
  } else if (str_detect(config[["suffix"]] , "_offset_add")) {
    offset_val <- as.numeric(str_match(config[["suffix"]] , "_offset_add(.+)")[2])
    correlation_pairs_NPP[ , correlation:=correlation+offset_val]
  }
  
  fwrite(correlation_pairs_NPP, file = output_correlation_file)
}
