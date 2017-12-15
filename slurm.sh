#!/bin/bash

#SBATCH -p priority
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -n 1
#SBATCH -e snakemake.err
#SBATCH -o snakemake.log
#SBATCH -J TSSseq-snakemake

snakemake -p -R `snakemake --list-code-changes` --latency-wait 300 --rerun-incomplete --cluster-config cluster.yaml --use-conda --jobs 999 --cluster "sbatch -p {cluster.queue} -n {cluster.n} -t {cluster.time} --mem-per-cpu={cluster.mem} -J {cluster.name} -e {cluster.err} -o {cluster.log}"
