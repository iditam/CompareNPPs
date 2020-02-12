#!/bin/bash

nohup snakemake -s snakemake_ES_with_random_zscores.smk --configfile CONFIG_NAME.yaml --rerun-incomplete  -j 500 --restart-times 1 --latency-wait 60 \
  --cluster-config snakemake_ES_zscore_icore.json --cluster "sbatch --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --time={cluster.time} --output={cluster.log} --error={cluster.log}" >> zscores_snakemake.out &
