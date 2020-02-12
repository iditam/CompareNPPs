#!/bin/bash

nohup snakemake -s snakemake_ES_with_random.smk -k --configfile CONFIG_NAME.yaml --rerun-incomplete \
  -j 500 --cluster-config snakemake_ES_with_random_icore.json --cluster "sbatch --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --time={cluster.time} --output={cluster.log} --error={cluster.log}" &
