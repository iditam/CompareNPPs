#!/bin/bash

nohup snakemake -s snakemake_correlations.smk --configfile config.yaml --rerun-incomplete -j 500 --latency-wait 60 --cluster-config snakemake_icore.json --cluster "sbatch --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --time={cluster.time} --output={cluster.log} --error={cluster.log}" &
