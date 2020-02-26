---
title: "Test2 - ES z-score for gene sets"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


**snakemake_ES_with_random.smk** - a snakemake pipeline that calculates ES for a single NPP, and the ES score for multiple random gene sets. The number of random gene sets is set by the "randomTimes" variable in the config file.

**run_snakemake_with_random.sh** - a script to run the pipeline "snakemake_ES_with_random.smk", specifying a config file and running the pipeline jobs on a cluster with SLURM job scheduling system. Replace "CONFIG_NAME.yaml" with correct config file name.
Run with:
> bash run_snakemake_with_random.sh

**snakemake_ES_with_random_zscores.smk** - a snakemake pipeline that includes the pipeline 'snakemake_ES_with_random.smk' and in addition steps for mirroring the random ES values and calculating ES z-score for each original gene set.

**run_snakemake_zscores.sh** - a script to run the pipeline "snakemake_ES_with_random_zscores.smk", specifying a config file and running the pipeline jobs on a cluster with SLURM job scheduling system. Replace "CONFIG_NAME.yaml" with correct config file name.
Run with:
> bash run_snakemake_zscores.sh

**config_XXX.yaml** - config files with definitions, used for the snakemake pipelines snakemake_ES_with_random.smk & snakemake_ES_with_random_zscores.smk

