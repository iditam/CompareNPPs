
**snakemake_correlations.smk** - a snakemake pipeline that creates a list with all gene pairs and the correlation between their profiles for each different table of NPP.

**run_snakemake_correlations.sh** - a script to run the pipeline "snakemake_correlations.smk", specifying a config file and running the pipeline jobs on a cluster with SLURM job scheduling system. Replace "CONFIG_NAME.yaml" with correct config file name.
Run with:
> bash run_snakemake_correlations.sh

**snakemake_ROC_multiple.smk** - a snakemake pipeline that includes the pipeline 'snakemake_correlations.smk' and in addition steps for calculating ROC curves for all NPP tables. The positive and negative datasets for the calculation are defined in the config file

**run_snakemake_ROC_multiple.sh** - a script to run the pipeline "snakemake_ROC_multiple.smk", specifying a config file and running the pipeline jobs on a cluster with SLURM job scheduling system. Replace "CONFIG_NAME.yaml" with correct config file name.
Run with:
> bash run_snakemake_ROC_multiple.sh

**config_XXX.yaml** - config files with definitions, used for the snakemake pipelines snakemake_correlations.smk & snakemake_ROC_multiple.smk

