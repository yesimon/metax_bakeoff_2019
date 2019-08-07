# metax

Metagenomic taxonomic classification benchmarking

The config.yaml file needs to be supplied but an example one is in `config/`. The ipython notebook
in `plotting/` is used to generate all plots and performance measurement. The `dropcache.c` file in
`scripts/` needs to be compiled and placed in the user's `$PATH` for rerun benchmarking. A set
conda package requirements are given by `conda-requirements.txt`.

Most of the workflow is implemented in Snakefile to let snakemake handle as much of the
benchmarking and dependencies as possible. Snakemake doesn't properly manage batching when sharing
a common resource easily, like a large metagenomic database wanting all of its dependent jobs run
serially so it won't have to reload the database from memory. To alleviate this and improve
performance, run separate executions of snakemake for each method, such as `snakemake
kraken_default_all; snakemake kraken_refseqc_all`.
